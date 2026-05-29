"""Centralized runtime configuration for genesis-bio-mcp.

All values can be overridden via environment variables with the ``GENESIS_``
prefix (e.g. ``GENESIS_HTTPX_TIMEOUT=60.0``) or via a ``.env`` file in the
working directory.

Usage::

    from genesis_bio_mcp.config.settings import settings

    client = httpx.AsyncClient(timeout=settings.httpx_timeout)
"""

from __future__ import annotations

from pathlib import Path

from pydantic import Field
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """Runtime settings for the genesis-bio-mcp server.

    All fields are read-only after construction. Change behaviour by setting
    the corresponding ``GENESIS_*`` environment variable before starting the
    server.

    Environment variables:
        GENESIS_HTTPX_TIMEOUT         — float, seconds (default 30.0)
        GENESIS_DEPMAP_CACHE_PATH     — path string (default data/depmap_cache.csv)
        GENESIS_DEPMAP_CACHE_MAX_AGE_DAYS — int (default 7)
        GENESIS_DEPMAP_TASK_TIMEOUT_SECS  — float (default 120.0)
        GENESIS_GWAS_CACHE_PATH       — path string (default data/gwas_cache.json)
        GENESIS_GWAS_CACHE_TTL_SECS   — int, seconds (default 86400 = 24 h)
        GENESIS_EFO_CACHE_PATH        — path string (default data/efo_cache.json)
        GENESIS_EFO_CACHE_TTL_SECS    — int, seconds (default 604800 = 7 days)
        GENESIS_CHEMBL_SEMAPHORE_LIMIT   — int (default 2)
        GENESIS_PUBCHEM_SEMAPHORE_LIMIT  — int (default 3)
        GENESIS_REACTOME_SEMAPHORE_LIMIT — int (default 3)
        GENESIS_STRING_REQUIRED_SCORE — int 0–1000 (default 700)
        GENESIS_STRING_NETWORK_LIMIT  — int (default 20)
        GENESIS_FOLDSEEK_DATABASE     — string (default afdb-swissprot)
        GENESIS_FOLDSEEK_POLL_TIMEOUT_SECS  — float (default 120.0)
        GENESIS_FOLDSEEK_POLL_INTERVAL_SECS — float (default 3.0)
        GENESIS_FOLDSEEK_MAX_HITS     — int (default 20)
        GENESIS_SABDAB_CACHE_PATH     — path string (default data/sabdab_cache.tsv)
        GENESIS_SABDAB_CACHE_TTL_SECS — int, seconds (default 604800 = 7 days)
        GENESIS_CLAUDE_MODEL          — string (default claude-sonnet-4-6)
    """

    model_config = SettingsConfigDict(
        env_prefix="GENESIS_",
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore",
    )

    # ---------------------------------------------------------------------------
    # HTTP client
    # ---------------------------------------------------------------------------

    httpx_timeout: float = Field(
        default=30.0,
        description="Default timeout in seconds for the shared httpx.AsyncClient.",
        gt=0,
    )

    # ---------------------------------------------------------------------------
    # DepMap cache
    # ---------------------------------------------------------------------------

    depmap_cache_path: Path = Field(
        default=Path("data/depmap_cache.csv"),
        description="Disk cache path for the DepMap gene dependency CSV.",
    )
    depmap_cache_max_age_days: int = Field(
        default=7,
        description="Maximum age in days before the DepMap disk cache is re-downloaded.",
        gt=0,
    )
    depmap_task_timeout_secs: float = Field(
        default=120.0,
        description="Timeout in seconds when polling a DepMap Celery async task.",
        gt=0,
    )

    # ---------------------------------------------------------------------------
    # GWAS cache
    # ---------------------------------------------------------------------------

    gwas_cache_path: Path = Field(
        default=Path("data/gwas_cache.json"),
        description="Disk cache path for GWAS Catalog association results.",
    )
    gwas_cache_ttl_secs: int = Field(
        default=86400,
        description="TTL in seconds for GWAS disk cache entries (default 24 h).",
        gt=0,
    )

    # ---------------------------------------------------------------------------
    # EFO ontology cache
    # ---------------------------------------------------------------------------

    efo_cache_path: Path = Field(
        default=Path("data/efo_cache.json"),
        description="Disk cache path for EFO OLS4 term resolution results.",
    )
    efo_cache_ttl_secs: int = Field(
        default=604800,
        description="TTL in seconds for EFO disk cache entries (default 7 days).",
        gt=0,
    )

    # ---------------------------------------------------------------------------
    # Semaphore limits — concurrent outbound API requests per service
    # ---------------------------------------------------------------------------

    chembl_semaphore_limit: int = Field(
        default=2,
        description="Max concurrent requests to ChEMBL (~1 req/s rate limit).",
        gt=0,
    )
    pubchem_semaphore_limit: int = Field(
        default=3,
        description="Max concurrent requests to PubChem.",
        gt=0,
    )
    reactome_semaphore_limit: int = Field(
        default=3,
        description="Max concurrent requests to Reactome AnalysisService / ContentService.",
        gt=0,
    )

    # ---------------------------------------------------------------------------
    # STRING protein-interaction network defaults
    # ---------------------------------------------------------------------------

    string_required_score: int = Field(
        default=700,
        description="Default STRING combined-score threshold (0–1000) for reported edges. "
        "700 = high confidence; lower to widen the network, raise to tighten it.",
        ge=0,
        le=1000,
    )
    string_network_limit: int = Field(
        default=20,
        description="Default maximum number of STRING interaction partners to return.",
        gt=0,
        le=200,
    )

    # ---------------------------------------------------------------------------
    # Foldseek structural-similarity search server
    # ---------------------------------------------------------------------------

    foldseek_database: str = Field(
        default="afdb-swissprot",
        description="Default Foldseek target database path (e.g. 'afdb-swissprot', 'pdb', "
        "'afdb-proteome'). See https://search.foldseek.com/api/databases.",
    )
    foldseek_poll_timeout_secs: float = Field(
        default=120.0,
        description="Max seconds to poll a Foldseek search ticket before giving up.",
        gt=0,
    )
    foldseek_poll_interval_secs: float = Field(
        default=3.0,
        description="Seconds between Foldseek ticket status polls.",
        gt=0,
    )
    foldseek_max_hits: int = Field(
        default=20,
        description="Default maximum number of structural-homolog hits to return.",
        gt=0,
        le=100,
    )

    # ---------------------------------------------------------------------------
    # SAbDab cache
    # ---------------------------------------------------------------------------

    sabdab_cache_path: Path = Field(
        default=Path("data/sabdab_cache.tsv"),
        description="Disk cache path for the SAbDab summary TSV.",
    )
    sabdab_cache_ttl_secs: int = Field(
        default=604800,
        description="TTL in seconds for the SAbDab disk cache (default 7 days).",
        gt=0,
    )

    # ---------------------------------------------------------------------------
    # GTEx tissue expression cache
    # ---------------------------------------------------------------------------

    gtex_cache_path: Path = Field(
        default=Path("data/gtex_cache.json"),
        description="Disk cache path for GTEx median gene-expression responses.",
    )
    gtex_cache_ttl_secs: int = Field(
        default=604800,
        description="TTL in seconds for the GTEx disk cache (default 7 days).",
        gt=0,
    )

    # ---------------------------------------------------------------------------
    # Human Protein Atlas cache
    # ---------------------------------------------------------------------------

    hpa_cache_path: Path = Field(
        default=Path("data/hpa_cache.json"),
        description="Disk cache path for HPA per-gene responses.",
    )
    hpa_cache_ttl_secs: int = Field(
        default=604800,
        description="TTL in seconds for the HPA disk cache (default 7 days).",
        gt=0,
    )

    # ---------------------------------------------------------------------------
    # OpenFDA drug-safety cache
    # ---------------------------------------------------------------------------

    openfda_cache_path: Path = Field(
        default=Path("data/openfda_cache.json"),
        description="Disk cache path for OpenFDA per-drug safety signals.",
    )
    openfda_cache_ttl_secs: int = Field(
        default=604800,
        description="TTL in seconds for the OpenFDA disk cache (default 7 days).",
        gt=0,
    )

    # ---------------------------------------------------------------------------
    # Workflow agent
    # ---------------------------------------------------------------------------

    claude_model: str = Field(
        default="claude-sonnet-4-6",
        description="Claude model ID used by run_biology_workflow.",
    )


#: Singleton settings instance. Import and use this directly.
settings = Settings()
