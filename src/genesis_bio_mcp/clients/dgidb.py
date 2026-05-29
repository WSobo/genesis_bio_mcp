"""DGIdb 5.0 GraphQL client for drug-gene interactions."""

from __future__ import annotations

import logging

import httpx

from genesis_bio_mcp.models import DrugInteraction

logger = logging.getLogger(__name__)

_DGIDB_URL = "https://dgidb.org/api/graphql"

# Interaction types that represent direct target engagement (inhibition, activation, binding).
# Other types (substrate, inducer, suppressor, etc.) are kept but sorted to the bottom.
_DIRECT_TYPES = {"inhibitor", "antagonist", "blocker", "agonist", "modulator", "binder"}

_QUERY = """
query GeneInteractions($gene: String!) {
  genes(names: [$gene]) {
    nodes {
      name
      interactions {
        drug {
          name
          approved
        }
        interactionScore
        interactionTypes {
          type
          directionality
        }
        interactionClaims {
          source {
            sourceDbName
          }
        }
      }
    }
  }
}
"""


class DGIdbClient:
    def __init__(self, client: httpx.AsyncClient) -> None:
        self._client = client

    async def get_drug_interactions(self, gene_symbol: str) -> list[DrugInteraction]:
        """Return drugs known to interact with the gene from DGIdb."""
        try:
            resp = await self._client.post(
                _DGIDB_URL,
                json={"query": _QUERY, "variables": {"gene": gene_symbol}},
                headers={"Content-Type": "application/json"},
                timeout=20.0,
            )
            resp.raise_for_status()
            body = resp.json()
            if "errors" in body:
                logger.warning("DGIdb GraphQL errors: %s", body["errors"])
                return []
            return _parse_interactions(body)
        except Exception as exc:
            logger.warning("DGIdb request failed for %s: %s", gene_symbol, exc)
            return []


def _max_opt(a: float | None, b: float | None) -> float | None:
    """Return the larger of two optional floats, ignoring ``None``."""
    vals = [v for v in (a, b) if v is not None]
    return max(vals) if vals else None


def _parse_interactions(body: dict) -> list[DrugInteraction]:
    nodes = body.get("data", {}).get("genes", {}).get("nodes", [])
    if not nodes:
        return []

    seen: dict[str, DrugInteraction] = {}
    for node in nodes:
        for interaction in node.get("interactions", []):
            drug = interaction.get("drug") or {}
            drug_name = (drug.get("name") or "").strip()
            if not drug_name:
                continue

            approved = bool(drug.get("approved", False))

            score_raw = interaction.get("interactionScore")
            score = float(score_raw) if isinstance(score_raw, (int, float)) else None

            # Determine interaction type from interactionTypes list
            interaction_types = interaction.get("interactionTypes") or []
            itype = None
            if interaction_types:
                itype = interaction_types[0].get("type")

            # Collect source names
            claims = interaction.get("interactionClaims") or []
            sources = list(
                {
                    (c.get("source") or {}).get("sourceDbName", "")
                    for c in claims
                    if (c.get("source") or {}).get("sourceDbName")
                }
            )

            # Estimate phase: approved → 4, else None (DGIdb doesn't always have phase)
            phase = 4 if approved else None

            key = drug_name.lower()
            if key not in seen:
                seen[key] = DrugInteraction(
                    drug_name=drug_name,
                    interaction_type=itype,
                    phase=phase,
                    approved=approved,
                    interaction_score=score,
                    sources=sorted(sources),
                )
            else:
                prev = seen[key]
                best_score = _max_opt(prev.interaction_score, score)
                if approved and not prev.approved:
                    # Upgrade to approved if we see a better record.
                    seen[key] = DrugInteraction(
                        drug_name=drug_name,
                        interaction_type=prev.interaction_type or itype,
                        phase=4,
                        approved=True,
                        interaction_score=best_score,
                        sources=sorted(set(prev.sources + sources)),
                    )
                elif best_score != prev.interaction_score:
                    # Same approval status, but keep the strongest score seen.
                    seen[key] = prev.model_copy(update={"interaction_score": best_score})

    # Sort: approved first, then direct interaction types, then alphabetical.
    # This surfaces confirmed inhibitors/modulators before substrate/inducer noise.
    ordered = sorted(
        seen.values(),
        key=lambda d: (
            not d.approved,
            d.interaction_type is None or d.interaction_type.lower() not in _DIRECT_TYPES,
            d.drug_name.lower(),
        ),
    )
    return _collapse_salt_forms(ordered)


def _collapse_salt_forms(drugs: list[DrugInteraction]) -> list[DrugInteraction]:
    """Merge records where one drug_name is a token-prefix of another.

    DGIdb reports salt forms as separate records (``"FILGOTINIB"`` and
    ``"FILGOTINIB MALEATE"``), double-counting the same INN. Pharma salt forms
    are universally formatted as ``"<INN> <counter-ion>"``, so a multi-token
    name whose first token matches a shorter single-token record's full name
    is virtually always the same molecule.

    Merges the longer (salt) record into the shorter (parent) record, unioning
    sources and preserving the strongest approved/phase signal. No hardcoded
    salt vocabulary — the decision is based purely on structural name prefix.
    """
    if len(drugs) < 2:
        return drugs

    by_name: dict[str, DrugInteraction] = {d.drug_name.lower(): d for d in drugs}

    # A record is a salt form of a "parent" when its first whitespace token
    # exactly matches a single-token record's full name. Restricting the parent
    # to single-token names avoids false positives on distinct multi-word drugs.
    merged_into: dict[str, str] = {}
    for name_l in list(by_name.keys()):
        tokens = name_l.split()
        if len(tokens) < 2:
            continue
        parent_key = tokens[0]
        if parent_key in by_name and parent_key != name_l:
            merged_into[name_l] = parent_key

    if not merged_into:
        return drugs

    for salt_key, parent_key in merged_into.items():
        salt = by_name[salt_key]
        parent = by_name[parent_key]
        by_name[parent_key] = DrugInteraction(
            drug_name=parent.drug_name,
            interaction_type=parent.interaction_type or salt.interaction_type,
            phase=max(parent.phase or 0, salt.phase or 0) or None,
            approved=parent.approved or salt.approved,
            interaction_score=_max_opt(parent.interaction_score, salt.interaction_score),
            sources=sorted(set(parent.sources) | set(salt.sources)),
        )

    return [by_name[d.drug_name.lower()] for d in drugs if d.drug_name.lower() not in merged_into]
