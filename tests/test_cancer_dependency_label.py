"""Regression test: the cancer-dependency table must not mislabel diseases as cell lines.

``get_cancer_dependency`` populates the ``cell_lines`` table from Open Targets
*disease-level* somatic-mutation evidence in every live path (the DepMap summary
endpoint exposes only counts, not per-line CERES). The model flags this with
``ceres_is_approximated=True``. The Markdown table header must reflect that — a
proxy table is a *disease* table, not a measured *cell line* table.
"""

from __future__ import annotations

from genesis_bio_mcp.models import CancerDependency, CellLineEssentiality


def _proxy_dependency(*, approximated: bool) -> CancerDependency:
    return CancerDependency(
        gene_symbol="BRAF",
        mean_ceres_score=-0.8,
        fraction_dependent_lines=0.4,
        pan_essential=False,
        top_dependent_lineages=["melanoma", "thyroid carcinoma"],
        cell_lines=[
            CellLineEssentiality(
                cell_line="melanoma",  # a DISEASE, not a cell line, in the proxy path
                lineage="skin cancer",
                ceres_score=-1.6,
                is_dependent=True,
            )
        ],
        data_source="Open Targets Platform v4 — somatic mutation evidence (proxy)",
        ceres_is_approximated=approximated,
    )


def test_proxy_table_labels_rows_as_diseases_not_cell_lines() -> None:
    md = _proxy_dependency(approximated=True).to_markdown()
    # The misleading column header is gone for proxy/disease-level data...
    assert "| Cell Line | Lineage |" not in md
    # ...and the honest disease/somatic-mutation framing is present.
    assert "Disease (somatic-mutation evidence)" in md
    # The disease name still renders in the row.
    assert "melanoma" in md


def test_real_per_line_table_keeps_cell_line_header() -> None:
    # If a genuine per-line CERES source is ever wired (not approximated),
    # the cell-line header is correct and must be preserved.
    md = _proxy_dependency(approximated=False).to_markdown()
    assert "| Cell Line | Lineage | CERES score | Dependent? |" in md
    assert "Disease (somatic-mutation evidence)" not in md
