"""ChEMBL bioactivity client.

Provides quantitative IC50/Ki/Kd potency data from the EMBL-EBI ChEMBL database.
No API key required; rate limit is ~1 req/sec (conservative semaphore at 2).

Flow per gene:
  1. Target search by gene symbol → extract first human single-protein target_chembl_id
  2. Activity query filtered to IC50/Ki/Kd with a pChEMBL value → top 20 by potency
"""

from __future__ import annotations

import asyncio
import logging

import httpx

from genesis_bio_mcp.config.settings import settings
from genesis_bio_mcp.models import ChEMBLActivity, ChEMBLCompounds, OffTarget

logger = logging.getLogger(__name__)

_BASE = "https://www.ebi.ac.uk/chembl/api/data"
_TARGET_SEARCH_URL = f"{_BASE}/target/search"
_MOLECULE_URL = f"{_BASE}/molecule"
_ACTIVITY_URL = f"{_BASE}/activity"
_ASSAY_URL = f"{_BASE}/assay"
_MECHANISM_URL = f"{_BASE}/mechanism"

_SEMAPHORE = asyncio.Semaphore(settings.chembl_semaphore_limit)


class ChEMBLClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client

    async def get_compounds(self, gene_symbol: str) -> ChEMBLCompounds | None:
        """Return quantitative bioactivity data for compounds active on gene_symbol's target."""
        target_id = await self._resolve_target(gene_symbol)
        if not target_id:
            logger.debug("ChEMBL: no human single-protein target found for %s", gene_symbol)
            return None

        # Activities (potency) and mechanisms (clinical precedent) are
        # independent endpoints — fetch them concurrently.
        activities, mechanism = await asyncio.gather(
            self._fetch_activities(target_id),
            self._fetch_mechanisms(target_id),
        )
        max_phase, moas, action_types = mechanism
        if not activities:
            return ChEMBLCompounds(
                gene_symbol=gene_symbol,
                target_chembl_id=target_id,
                total_active_compounds=0,
                best_pchembl=None,
                max_clinical_phase=max_phase,
                mechanisms_of_action=moas,
                action_types=action_types,
                compounds=[],
            )

        best = max(a.pchembl_value for a in activities)
        # Split bests by assay type so the scoring axis can prefer functional/
        # cell-based potency over binding-only (which over-credits scaffold
        # binders for undruggable targets like MYC).
        functional = [a.pchembl_value for a in activities if (a.assay_type or "").upper() == "F"]
        binding = [a.pchembl_value for a in activities if (a.assay_type or "").upper() == "B"]
        return ChEMBLCompounds(
            gene_symbol=gene_symbol,
            target_chembl_id=target_id,
            total_active_compounds=len(activities),
            best_pchembl=round(best, 2),
            best_pchembl_functional=round(max(functional), 2) if functional else None,
            best_pchembl_binding=round(max(binding), 2) if binding else None,
            max_clinical_phase=max_phase,
            mechanisms_of_action=moas,
            action_types=action_types,
            compounds=activities[:20],
        )

    async def _resolve_target(self, gene_symbol: str) -> str | None:
        """Find the ChEMBL target ID for a human gene symbol."""
        try:
            async with _SEMAPHORE:
                resp = await self._client.get(
                    _TARGET_SEARCH_URL,
                    params={"q": gene_symbol, "format": "json", "limit": 10},
                    timeout=20.0,
                )
                resp.raise_for_status()

            data = resp.json()
            for target in data.get("targets") or []:
                # Filter: single protein, Homo sapiens
                if (
                    target.get("target_type") == "SINGLE PROTEIN"
                    and target.get("organism") == "Homo sapiens"
                ):
                    return target.get("target_chembl_id")
            return None

        except Exception as exc:
            logger.warning("ChEMBL target search failed for %s: %s", gene_symbol, exc)
            return None

    async def _fetch_activities(self, target_id: str) -> list[ChEMBLActivity]:
        """Fetch IC50/Ki/Kd activities with a pChEMBL value, sorted by potency."""
        try:
            async with _SEMAPHORE:
                resp = await self._client.get(
                    _ACTIVITY_URL,
                    params={
                        "target_chembl_id": target_id,
                        "standard_type__in": "IC50,Ki,Kd,EC50",
                        "pchembl_value__isnull": "false",
                        "limit": 100,
                        "format": "json",
                    },
                    timeout=30.0,
                )
                resp.raise_for_status()

            data = resp.json()
            raw = data.get("activities") or []

            activities: list[ChEMBLActivity] = []
            seen: set[str] = set()
            assay_by_mol: dict[
                str, str
            ] = {}  # molecule → assay_chembl_id (for confidence enrichment)

            for a in raw:
                mol_id = a.get("molecule_chembl_id", "")
                pchembl = a.get("pchembl_value")
                if not mol_id or pchembl is None:
                    continue
                try:
                    pchembl_f = float(pchembl)
                except (ValueError, TypeError):
                    continue

                # Deduplicate by molecule: keep best pChEMBL value per compound
                if mol_id in seen:
                    continue
                seen.add(mol_id)

                assay_id = a.get("assay_chembl_id")
                if assay_id:
                    assay_by_mol[mol_id] = assay_id

                activities.append(
                    ChEMBLActivity(
                        molecule_chembl_id=mol_id,
                        molecule_name=a.get("molecule_pref_name"),
                        standard_type=a.get("standard_type", ""),
                        pchembl_value=round(pchembl_f, 2),
                        assay_description=a.get("assay_description", "")[:120]
                        if a.get("assay_description")
                        else None,
                        assay_type=a.get("assay_type") or None,
                        # ChEMBL's /activity rows expose the species via
                        # ``target_organism``; ``assay_organism`` is almost
                        # always null. Fall back to the latter for assays where
                        # the assay system differs from the target species
                        # (e.g. human protein expressed in insect cells).
                        assay_organism=a.get("target_organism") or a.get("assay_organism") or None,
                        assay_cell_type=a.get("assay_cell_type") or None,
                        bao_format=a.get("bao_label") or a.get("bao_format") or None,
                        # /activity carries NO confidence_score (it's on /assay) — enriched below.
                        confidence_score=None,
                    )
                )

            # ChEMBL's /activity endpoint has no confidence_score field — it lives on /assay.
            # Enrich so the low-target-confidence caveat in to_markdown actually fires.
            conf_map = await self._fetch_assay_confidence(set(assay_by_mol.values()))
            for act in activities:
                aid = assay_by_mol.get(act.molecule_chembl_id)
                if aid and aid in conf_map:
                    act.confidence_score = conf_map[aid]

            # Sort by potency descending
            activities.sort(key=lambda x: x.pchembl_value, reverse=True)
            return activities

        except Exception as exc:
            logger.warning("ChEMBL activity fetch failed for %s: %s", target_id, exc)
            return []

    async def _fetch_assay_confidence(self, assay_ids: set[str]) -> dict[str, int]:
        """Map assay_chembl_id → confidence_score from /assay (batched). Never raises."""
        ids = sorted(a for a in assay_ids if a)
        out: dict[str, int] = {}
        for i in range(0, len(ids), 50):
            chunk = ids[i : i + 50]
            try:
                async with _SEMAPHORE:
                    resp = await self._client.get(
                        _ASSAY_URL,
                        params={
                            "assay_chembl_id__in": ",".join(chunk),
                            "format": "json",
                            "limit": len(chunk),
                        },
                        timeout=20.0,
                    )
                    resp.raise_for_status()
            except Exception as exc:
                logger.warning("ChEMBL assay confidence fetch failed: %s", exc)
                continue
            for assay in resp.json().get("assays") or []:
                cs = assay.get("confidence_score")
                aid = assay.get("assay_chembl_id")
                if aid and cs is not None:
                    try:
                        out[aid] = int(cs)
                    except (TypeError, ValueError):
                        pass
        return out

    async def _fetch_mechanisms(self, target_id: str) -> tuple[float | None, list[str], list[str]]:
        """Fetch mechanism-of-action records for the target.

        Returns ``(max_clinical_phase, mechanisms_of_action, action_types)``.
        The ChEMBL ``/mechanism`` endpoint carries ``max_phase`` inline, so the
        most advanced clinical stage of any drug acting on this target comes
        from a single call (no per-molecule fan-out). Never raises.
        """
        try:
            async with _SEMAPHORE:
                resp = await self._client.get(
                    _MECHANISM_URL,
                    params={
                        "target_chembl_id": target_id,
                        "limit": 100,
                        "format": "json",
                    },
                    timeout=20.0,
                )
                resp.raise_for_status()
            rows = resp.json().get("mechanisms") or []
        except Exception as exc:
            logger.warning("ChEMBL mechanism fetch failed for %s: %s", target_id, exc)
            return None, [], []

        max_phase: float | None = None
        moas: list[str] = []
        action_types: list[str] = []
        for m in rows:
            phase_raw = m.get("max_phase")
            try:
                phase = float(phase_raw) if phase_raw is not None else None
            except (TypeError, ValueError):
                phase = None
            if phase is not None and (max_phase is None or phase > max_phase):
                max_phase = phase
            moa = m.get("mechanism_of_action")
            if isinstance(moa, str) and moa and moa not in moas:
                moas.append(moa)
            at = m.get("action_type")
            if isinstance(at, str) and at and at not in action_types:
                action_types.append(at)
        return max_phase, moas, action_types

    async def _get_retry(
        self, url: str, *, params: dict | None = None, timeout: float = 30.0, attempts: int = 3
    ) -> httpx.Response:
        """GET with backoff retry on transient ChEMBL 5xx / network errors (4xx returned as-is).

        ChEMBL intermittently 500s / times out under load. Retrying transient failures keeps the
        off-target lookup robust. Raises the last error after ``attempts`` tries.
        """
        last: Exception | None = None
        for i in range(attempts):
            try:
                async with _SEMAPHORE:
                    resp = await self._client.get(url, params=params, timeout=timeout)
                if resp.status_code < 500:
                    return resp
                last = httpx.HTTPStatusError(
                    f"server {resp.status_code}", request=resp.request, response=resp
                )
            except Exception as exc:
                last = exc
            if i < attempts - 1:
                await asyncio.sleep(1.5 * (i + 1))
        raise last  # type: ignore[misc]

    async def get_off_targets(self, inchikey: str) -> tuple[str | None, list[OffTarget]]:
        """Measured off-target profile: every human target a compound has pChEMBL activity against.

        ``inchikey`` is the standard InChIKey of the (standardized) query molecule. Resolves it to a
        ChEMBL molecule and aggregates all its human bioactivities by target (using the activity
        rows' ``target_pref_name`` — ChEMBL's ``/target`` endpoint is too unreliable to enrich gene
        symbols here; the corpus predicted layer supplies those). Never raises — returns ``(None, [])``
        if the structure is not in ChEMBL or on any failure. Off-targets sorted by best pChEMBL desc.
        """
        mol_id = await self._resolve_molecule_id(inchikey)
        if not mol_id:
            return None, []
        agg = await self._fetch_molecule_activities(mol_id)
        off = [
            OffTarget(
                target_chembl_id=tid,
                target_pref_name=a["pref_name"],
                best_pchembl=round(a["best"], 2),
                n_activities=a["n"],
            )
            for tid, a in agg.items()
        ]
        off.sort(key=lambda o: o.best_pchembl, reverse=True)
        return mol_id, off

    async def _resolve_molecule_id(self, inchikey: str) -> str | None:
        """Resolve a standard InChIKey to a ChEMBL molecule ID via the resource endpoint.

        Uses ``/molecule/{inchikey}.json`` (the resource path) — the
        ``?molecule_structures__standard_inchi_key=`` *filter* 500s server-side. Never raises.
        """
        try:
            resp = await self._get_retry(f"{_MOLECULE_URL}/{inchikey}.json", timeout=20.0)
            if resp.status_code == 404:
                return None
            resp.raise_for_status()
            return resp.json().get("molecule_chembl_id")
        except Exception as exc:
            logger.warning("ChEMBL molecule lookup failed for %s: %s", inchikey, exc)
            return None

    async def _fetch_molecule_activities(self, mol_id: str) -> dict[str, dict]:
        """Aggregate a molecule's human pChEMBL activities by target. Never raises.

        Returns ``{target_chembl_id: {"pref_name", "best", "n"}}``. Pages by explicit offset
        (page size 200 — larger limits 500 under load; cap 10 pages) to capture a promiscuous
        compound's full target set without next-URL parsing.
        """
        page_size = 200
        agg: dict[str, dict] = {}
        try:
            for page in range(10):
                resp = await self._get_retry(
                    _ACTIVITY_URL,
                    params={
                        "molecule_chembl_id": mol_id,
                        "pchembl_value__isnull": "false",
                        "limit": page_size,
                        "offset": page * page_size,
                        "format": "json",
                    },
                    timeout=30.0,
                )
                resp.raise_for_status()
                data = resp.json()
                rows = data.get("activities") or []
                for a in rows:
                    if a.get("target_organism") != "Homo sapiens":
                        continue
                    tid = a.get("target_chembl_id")
                    pchembl = a.get("pchembl_value")
                    if not tid or pchembl is None:
                        continue
                    try:
                        pf = float(pchembl)
                    except (TypeError, ValueError):
                        continue
                    e = agg.setdefault(
                        tid, {"pref_name": a.get("target_pref_name") or tid, "best": pf, "n": 0}
                    )
                    e["best"] = max(e["best"], pf)
                    e["n"] += 1
                total = (data.get("page_meta") or {}).get("total_count") or 0
                if not rows or (page + 1) * page_size >= total:
                    break
        except Exception as exc:
            logger.warning("ChEMBL molecule-activity fetch failed for %s: %s", mol_id, exc)
        return agg
