"""Shared mock fixtures for genesis-bio-mcp tests."""

import httpx
import pytest

from genesis_bio_mcp.models import (
    CancerDependency,
    CellLineEssentiality,
    CompoundActivity,
    Compounds,
    DiseaseLinkEvidence,
    GeneResolution,
    GwasEvidence,
    GwasHit,
    KnownVariant,
    ProteinInfo,
    TargetDiseaseAssociation,
)

# ---------------------------------------------------------------------------
# Shared HTTP client fixture
# ---------------------------------------------------------------------------


@pytest.fixture
def http_client():
    return httpx.AsyncClient()


# ---------------------------------------------------------------------------
# Raw API response mocks (match real API shapes)
# ---------------------------------------------------------------------------


MOCK_UNIPROT_BRAF = {
    "primaryAccession": "P15056",
    "entryType": "UniProtKB reviewed (Swiss-Prot)",
    "genes": [
        {
            "geneName": {"value": "BRAF"},
            "synonyms": [{"value": "BRAF1"}, {"value": "RAFB1"}],
        }
    ],
    "proteinDescription": {
        "recommendedName": {"fullName": {"value": "Serine/threonine-protein kinase B-raf"}}
    },
    "organism": {"scientificName": "Homo sapiens"},
    "comments": [
        {
            "commentType": "FUNCTION",
            "texts": [
                {
                    "value": (
                        "Protein kinase involved in the transduction of mitogenic signals "
                        "from the cell membrane to the nucleus."
                    )
                }
            ],
        },
        {
            "commentType": "SUBCELLULAR LOCATION",
            "subcellularLocations": [
                {"location": {"value": "Cytoplasm"}},
                {"location": {"value": "Cell membrane"}},
            ],
        },
        {
            "commentType": "DISEASE",
            "disease": {"diseaseName": "Melanoma"},
        },
    ],
    "uniProtKBCrossReferences": [
        {"database": "PDB", "id": "1UWH"},
        {"database": "PDB", "id": "4MNE"},
        {
            "database": "Reactome",
            "id": "R-HSA-5673001",
            "properties": [{"key": "PathwayName", "value": "RAF/MAP kinase cascade"}],
        },
    ],
    "features": [
        {
            "type": "Natural variant",
            "location": {"start": {"value": 600}, "end": {"value": 600}},
            "alternativeSequence": {
                "originalSequence": "V",
                "alternativeSequences": ["E"],
            },
            "description": "In melanoma; somatic mutation",
        },
        {
            "type": "Disulfide bond",
            "location": {"start": {"value": 157}, "end": {"value": 162}},
            "description": "",
        },
    ],
}

# Real UniProt FASTA format — single-line header with OS/OX/GN metadata,
# sequence wrapped to 60 characters. Truncated to a manageable test size.
MOCK_UNIPROT_FASTA_BRAF = (
    ">sp|P15056|BRAF_HUMAN Serine/threonine-protein kinase B-raf OS=Homo sapiens OX=9606 GN=BRAF PE=1 SV=4\n"
    "MAALSGGGGGGAEPGQALFNGDMEPEAGAGAGAAASSAADPAIPEEVWNIKQMIKLTQEHI\n"
    "EALLDKFGGEHNPPSIYLEAYEEYTSKLDALQQREQQLLESLGNGTDFSVSSSASMDTVTS\n"
    "SSSSSLSVLPSSLSVFQNPTDVARSNPKSPQKPIVRVFLPNKQRTVVPARCGVTVRDSLKK\n"
)

MOCK_UNIPROT_FASTA_BRAF_SEQUENCE = (
    "MAALSGGGGGGAEPGQALFNGDMEPEAGAGAGAAASSAADPAIPEEVWNIKQMIKLTQEHI"
    "EALLDKFGGEHNPPSIYLEAYEEYTSKLDALQQREQQLLESLGNGTDFSVSSSASMDTVTS"
    "SSSSSLSVLPSSLSVFQNPTDVARSNPKSPQKPIVRVFLPNKQRTVVPARCGVTVRDSLKK"
)

MOCK_OT_GENE_SEARCH = {
    "data": {"search": {"hits": [{"id": "ENSG00000157764", "name": "BRAF", "entity": "target"}]}}
}

MOCK_OT_DISEASE_SEARCH = {
    "data": {
        "search": {
            "hits": [{"id": "EFO_0000389", "name": "cutaneous melanoma", "entity": "disease"}]
        }
    }
}

MOCK_OT_ASSOCIATION = {
    "data": {
        "target": {
            "associatedDiseases": {
                "count": 1,
                "rows": [
                    {
                        "score": 0.89,
                        "datatypeScores": [
                            {"id": "genetic_association", "score": 0.72},
                            {"id": "somatic_mutation", "score": 0.95},
                            {"id": "known_drug", "score": 0.88},
                            {"id": "literature", "score": 0.61},
                        ],
                        "datasourceScores": [],
                        "disease": {"id": "EFO_0000389", "name": "cutaneous melanoma"},
                    }
                ],
            }
        }
    }
}

# DepMap now uses Open Targets cancer association data as a proxy
MOCK_DEPMAP_OT_CANCER = {
    "data": {
        "target": {
            "associatedDiseases": {
                "count": 3,
                "rows": [
                    {
                        "score": 0.82,
                        "datatypeScores": [
                            {"id": "somatic_mutation", "score": 0.80},
                            {"id": "known_drug", "score": 0.88},
                        ],
                        "disease": {
                            "id": "EFO_0000756",
                            "name": "melanoma",
                            "therapeuticAreas": [{"name": "cancer or benign tumor"}],
                        },
                    },
                    {
                        "score": 0.75,
                        "datatypeScores": [
                            {"id": "somatic_mutation", "score": 0.80},
                        ],
                        "disease": {
                            "id": "EFO_0001253",
                            "name": "thyroid carcinoma",
                            "therapeuticAreas": [{"name": "cancer or benign tumor"}],
                        },
                    },
                    {
                        "score": 0.30,
                        "datatypeScores": [
                            {"id": "somatic_mutation", "score": 0.10},
                        ],
                        "disease": {
                            "id": "EFO_0000270",
                            "name": "asthma",
                            "therapeuticAreas": [{"name": "respiratory disease"}],
                        },
                    },
                ],
            }
        }
    }
}

MOCK_GWAS_SNP_RESPONSE = {
    "_embedded": {
        "singleNucleotidePolymorphisms": []  # Empty to avoid follow-up requests in unit tests
    }
}

MOCK_GWAS_ASSOCIATION_RESPONSE = {
    "_embedded": {
        "associations": [
            {
                "pvalue": 3.2e-15,
                "betaNum": 0.43,
                "efoTraits": [{"trait": "melanoma"}],
                "loci": [
                    {
                        "strongestRiskAlleles": [{"riskAlleleName": "rs1801516-A"}],
                        "authorReportedGenes": [{"geneName": "BRAF"}],
                    }
                ],
                "study": {
                    "studyAccession": "GCST001234",
                    "pubmedId": "12345678",
                    "initialSampleSize": "10000 European ancestry cases",
                },
                "_links": {
                    "study": {"href": "https://www.ebi.ac.uk/gwas/rest/api/studies/GCST001234"}
                },
            }
        ]
    }
}

MOCK_PUBCHEM_AIDS = {"IdentifierList": {"AID": [1259398, 686978, 504329]}}

MOCK_PUBCHEM_BIOACTIVITIES = {
    "Table": {
        "Columns": {"Column": ["CID", "Compound", "Activity Outcome", "Activity Value [uM]"]},
        "Row": [
            {"Cell": ["44462760", "Vemurafenib", "Active", "0.031"]},
            {"Cell": ["11338033", "Dabrafenib", "Active", "0.006"]},
            {"Cell": ["99999999", "Inactive Cpd", "Inactive", ""]},
        ],
    }
}

# New PubChem mocks for updated client (Entrez → active CIDs → properties)
MOCK_ENTREZ_AIDS = {
    "esearchresult": {
        "idlist": ["1259398", "686978"],
        "count": "2",
    }
}

MOCK_PUBCHEM_ACTIVE_CIDS = {"IdentifierList": {"CID": [44462760, 11338033]}}

# PUG REST primary path: assay/target/genesymbol/{symbol}/aids/JSON
MOCK_PUBCHEM_GENESYMBOL_AIDS = {"IdentifierList": {"AID": [1259398, 686978]}}

MOCK_PUBCHEM_PROPERTIES = {
    "PropertyTable": {
        "Properties": [
            {
                "CID": 44462760,
                "MolecularFormula": "C23H18ClF2N3O3S",
                "MolecularWeight": 489.9,
                "IUPACName": "vemurafenib",
            },
            {
                "CID": 11338033,
                "MolecularFormula": "C23H20F3N5O2S2",
                "MolecularWeight": 519.6,
                "IUPACName": "dabrafenib",
            },
        ]
    }
}


# ---------------------------------------------------------------------------
# Model-level fixtures (for e2e tests)
# ---------------------------------------------------------------------------


def build_mock_resolution(symbol: str = "BRAF") -> GeneResolution:
    return GeneResolution(
        hgnc_symbol=symbol,
        ncbi_gene_id="673",
        uniprot_accession="P15056",
        synonyms=["BRAF1", "RAFB1"],
        source="uniprot",
    )


def build_mock_protein_info(symbol: str = "BRAF") -> ProteinInfo:
    return ProteinInfo(
        uniprot_accession="P15056",
        gene_symbol=symbol,
        protein_name="Serine/threonine-protein kinase B-raf",
        organism="Homo sapiens",
        function_summary="Protein kinase involved in the transduction of mitogenic signals.",
        subcellular_locations=["Cytoplasm", "Cell membrane"],
        pathways=["RAF/MAP kinase cascade", "MAPK signaling pathway"],
        disease_associations=["Melanoma", "Lung adenocarcinoma"],
        pdb_structures=["1UWH", "4MNE"],
        known_variants=[
            KnownVariant(position="600", original="V", variant="E", disease="Melanoma")
        ],
        reviewed=True,
    )


def build_mock_association(
    gene: str = "BRAF", disease: str = "melanoma", score: float = 0.89
) -> TargetDiseaseAssociation:
    return TargetDiseaseAssociation(
        gene_symbol=gene,
        disease_name=disease,
        disease_efo_id="EFO_0000389",
        ensembl_id="ENSG00000157764",
        overall_score=score,
        genetic_association_score=0.72,
        somatic_mutation_score=0.95,
        known_drug_score=0.88,
        literature_mining_score=0.61,
        evidence_count=312,
        evidence_breakdown=[
            DiseaseLinkEvidence(evidence_type="somatic_mutation", score=0.95),
            DiseaseLinkEvidence(evidence_type="known_drug", score=0.88),
        ],
    )


def build_mock_dependency(
    gene: str = "BRAF",
    mean_score: float = -0.72,
    fraction_dependent: float = 0.61,
) -> CancerDependency:
    return CancerDependency(
        gene_symbol=gene,
        mean_ceres_score=mean_score,
        fraction_dependent_lines=fraction_dependent,
        pan_essential=False,
        top_dependent_lineages=["Skin", "Thyroid", "Colon"],
        cell_lines=[
            CellLineEssentiality(
                cell_line="A375", lineage="Skin", ceres_score=-1.82, is_dependent=True
            ),
            CellLineEssentiality(
                cell_line="SKMEL28",
                lineage="Skin",
                ceres_score=-1.65,
                is_dependent=True,
            ),
        ],
        data_source="DepMap Public 24Q4",
    )


def build_mock_gwas(gene: str = "BRAF", trait: str = "melanoma", n_hits: int = 3) -> GwasEvidence:
    hits = [
        GwasHit(
            study_accession=f"GCST{1000 + i}",
            trait=trait,
            mapped_gene=gene,
            risk_allele=f"rs{1000000 + i}-A",
            p_value=10 ** -(15 + i),
            beta_or_or=0.4,
            pubmed_id=f"1234567{i}",
        )
        for i in range(n_hits)
    ]
    return GwasEvidence(
        gene_symbol=gene,
        trait_query=trait,
        total_associations=n_hits,
        associations=hits,
        strongest_p_value=hits[0].p_value if hits else None,
    )


def build_mock_compounds(gene: str = "BRAF", n_active: int = 312) -> Compounds:
    return Compounds(
        gene_symbol=gene,
        total_active_compounds=n_active,
        compounds=[
            CompoundActivity(
                cid=44462760,
                name="Vemurafenib",
                molecular_formula="C23H18ClF2N3O3S",
                molecular_weight=489.9,
                activity_outcome="Active",
                activity_value=31.0,
                activity_type="IC50",
                assay_id=1259398,
            ),
            CompoundActivity(
                cid=11338033,
                name="Dabrafenib",
                molecular_formula="C23H20F3N5O2S2",
                molecular_weight=519.6,
                activity_outcome="Active",
                activity_value=6.0,
                activity_type="IC50",
                assay_id=1259398,
            ),
        ],
    )
