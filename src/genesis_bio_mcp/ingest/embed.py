"""ESM-2 protein-sequence embeddings (CPU, offline ingestion only).

``torch`` and ``fair-esm`` are OPTIONAL — installed via ``uv sync --group ingest`` — and
imported lazily so the serving install and the ``corpus_*`` tools never load an ML model.
We use ESM-2-150M (``esm2_t30_150M_UR50D``, 640-dim) mean-pooled: small enough to run on
CPU in a reasonable one-time batch, big enough to cluster a protein family meaningfully.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np

logger = logging.getLogger(__name__)

ESM2_MODEL = "esm2_t30_150M_UR50D"
EMBED_DIM = 640
_REPR_LAYER = 30  # final representation layer for the t30 (150M) model


def _require_esm():
    """Import torch + esm lazily, with an actionable error if the ingest extras are missing."""
    try:
        import esm
        import torch
    except ImportError as exc:  # pragma: no cover - exercised only without the extras
        raise ImportError(
            "ESM-2 embedding needs the optional 'ingest' dependencies (torch, fair-esm). "
            "Install them with:  uv sync --group ingest"
        ) from exc
    return torch, esm


def _mean_pool(token_representations: np.ndarray, seq_len: int) -> list[float]:
    """Mean-pool per-residue representations over the real residues only.

    fair-esm prepends a BOS token (index 0) and appends EOS (index ``seq_len + 1``), so the
    residue representations are rows ``1 .. seq_len``. Returns a length-``EMBED_DIM`` list.
    """
    import numpy as np

    arr = np.asarray(token_representations, dtype="float32")
    residues = arr[1 : seq_len + 1]
    if residues.shape[0] == 0:  # degenerate (empty/odd input) — fall back to all tokens
        residues = arr
    return residues.mean(axis=0).tolist()


def embed_sequences(sequences: list[str], *, batch_size: int = 8) -> list[list[float]]:
    """Compute mean-pooled ESM-2-150M embeddings for ``sequences`` on CPU.

    Returns one length-640 float list per input sequence, in order. Downloads the model
    weights (~600 MB) from a public URL on first use — no API key. Raises ``ImportError`` if
    the optional ``ingest`` extras are not installed.
    """
    torch, esm = _require_esm()
    model, alphabet = esm.pretrained.esm2_t30_150M_UR50D()
    model.eval()
    batch_converter = alphabet.get_batch_converter()

    out: list[list[float]] = []
    with torch.no_grad():
        for start in range(0, len(sequences), batch_size):
            chunk = sequences[start : start + batch_size]
            _, _, tokens = batch_converter([(str(i), s) for i, s in enumerate(chunk)])
            reps = model(tokens, repr_layers=[_REPR_LAYER])["representations"][_REPR_LAYER]
            for j, seq in enumerate(chunk):
                out.append(_mean_pool(reps[j].cpu().numpy(), len(seq)))
            logger.info(
                "Embedded %d/%d sequences", min(start + batch_size, len(sequences)), len(sequences)
            )
    return out
