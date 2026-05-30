-- genesis-bio-mcp v0.6.0 hybrid-retrieval corpus schema (Postgres + pgvector).
--
-- Idempotent: safe to run on every server start (CREATE ... IF NOT EXISTS).
-- Scope = one target family (the human kinome). Embedding columns are populated
-- by the OFFLINE ingestion tier; the serving tier only reads them.
--
-- Vector indexes are deliberately NOT created here: at corpus scale (~500 targets,
-- tens of thousands of compounds) exact kNN is sub-millisecond and avoids the
-- filtered-HNSW recall trap. Add HNSW only if the corpus grows (see docs/ROADMAP.md).

CREATE EXTENSION IF NOT EXISTS vector;

-- Single-row build manifest (FAIR provenance: what was indexed, from which releases).
CREATE TABLE IF NOT EXISTS corpus_manifest (
    id                       integer PRIMARY KEY DEFAULT 1,
    built_at                 timestamptz NOT NULL,
    target_family            text NOT NULL,
    chembl_release           text,
    uniprot_snapshot         text,
    protein_embedding_model  text,
    chem_embedding_model     text,
    target_count             integer NOT NULL DEFAULT 0,
    compound_count           integer NOT NULL DEFAULT 0,
    activity_count           integer NOT NULL DEFAULT 0,
    CONSTRAINT corpus_manifest_singleton CHECK (id = 1)
);

-- Keyed on UniProt accession — the FAIR cross-source identifier. ChEMBL's own
-- target id is kept as a (nullable) column and filled when activities are ingested.
CREATE TABLE IF NOT EXISTS targets (
    uniprot_accession   text PRIMARY KEY,
    target_chembl_id    text,
    gene_symbol         text,
    pref_name           text,
    organism            text,
    sequence            text,
    kinase_group        text,
    sequence_embedding  vector(640),    -- ESM-2-150M (esm2_t30_150M_UR50D), mean-pooled (offline, CPU)
    source              text,
    source_version      text,
    retrieved_at        timestamptz
);

CREATE TABLE IF NOT EXISTS compounds (
    molecule_chembl_id  text PRIMARY KEY,
    canonical_smiles    text,
    inchi               text,
    inchikey            text,
    mol_weight          double precision,
    morgan_fp           bit(2048),       -- ECFP4 (RDKit) — PRIMARY chemical similarity
    chem_embedding      vector(768),     -- ChemBERTa (optional/stretch)
    source              text,
    source_version      text,
    retrieved_at        timestamptz
);

CREATE TABLE IF NOT EXISTS activities (
    activity_id            bigint PRIMARY KEY,
    molecule_chembl_id     text REFERENCES compounds(molecule_chembl_id),
    uniprot_accession      text REFERENCES targets(uniprot_accession),
    target_chembl_id       text,          -- ChEMBL's own target id (provenance)
    standard_type          text,          -- IC50 / Ki / Kd
    standard_value         double precision,
    standard_units         text,
    pchembl_value          double precision,
    assay_chembl_id        text,
    assay_confidence_score integer,
    doc_chembl_id          text,          -- provenance: source publication
    source                 text,
    source_version         text,
    retrieved_at           timestamptz
);

CREATE INDEX IF NOT EXISTS idx_activities_target   ON activities (uniprot_accession);
CREATE INDEX IF NOT EXISTS idx_activities_molecule ON activities (molecule_chembl_id);
CREATE INDEX IF NOT EXISTS idx_activities_pchembl  ON activities (pchembl_value);
CREATE INDEX IF NOT EXISTS idx_targets_gene        ON targets (gene_symbol);
CREATE INDEX IF NOT EXISTS idx_compounds_inchikey  ON compounds (inchikey);
