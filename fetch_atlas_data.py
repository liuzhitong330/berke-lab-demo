"""
Fetch Allen Brain Atlas ISH gene expression data for dopamine pathway genes
across mouse striatal subregions, to illustrate the molecular substrate of
the dorsoventral reward-gradient described in Berke lab's work.

Data source: Allen Mouse Brain Atlas (2004), via public API at api.brain-map.org
All data is publicly available, no authentication required.
"""
import requests, json, time, csv
from pathlib import Path

OUT = Path(__file__).parent

# ── Striatal structures of interest ──────────────────────────────────────────
# Dorsal → Ventral ordering maps to extended → immediate reward time horizons
# (the axis studied in Berke lab's striatal gradient papers)
STRUCTURE_ACRONYMS = ["CP", "ACB", "OT"]
STRUCTURE_LABELS   = {
    "CP":  "Caudoputamen\n(dorsal striatum)",
    "ACB": "Nucleus accumbens\n(ventral striatum)",
    "OT":  "Olfactory tubercle\n(most ventral)"
}

# Dopamine pathway gene panel
GENES = [
    "Drd1",    # D1 receptor – direct pathway, reward
    "Drd2",    # D2 receptor – indirect pathway, aversion
    "Drd3",    # D3 receptor – mainly ventral
    "Slc6a3",  # DAT – dopamine reuptake transporter
    "Slc18a2", # VMAT2 – vesicular monoamine transporter
    "Th",      # Tyrosine hydroxylase – dopamine synthesis
    "Ppp1r1b", # DARPP-32 – dopaminergic signalling integrator
    "Penk",    # Enkephalin – indirect pathway marker
    "Pdyn",    # Dynorphin – direct pathway marker
    "Adora2a", # Adenosine A2a receptor – indirect pathway co-marker
]

BASE = "http://api.brain-map.org/api/v2/data/query.json"


def get_structure_ids(acronyms):
    """Return {acronym: id} for given Allen structure acronyms (mouse, ontology 1)."""
    acr_str = ",".join(f"'{a}'" for a in acronyms)
    url = (f"{BASE}?criteria=model::Structure,rma::criteria,"
           f"[acronym$in{acr_str}],"
           f"ontology[id$eq1]"
           f"&include=ontology"
           f"&num_rows=50")
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    data = r.json()["msg"]
    return {s["acronym"]: s["id"] for s in data}


def get_ish_datasets(gene_acronym, product_id=1):
    """Return list of SectionDataSet IDs for a gene in the Allen Mouse Brain Atlas."""
    url = (f"{BASE}?criteria=model::SectionDataSet,rma::criteria,"
           f"products[id$eq{product_id}],"
           f"genes[acronym$eq'{gene_acronym}']"
           f"&num_rows=5")
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return [d["id"] for d in r.json()["msg"]]


def get_structure_expression(dataset_id, structure_ids):
    """Return expression_energy per structure for one ISH dataset."""
    id_str = ",".join(str(i) for i in structure_ids.values())
    url = (f"{BASE}?criteria=model::StructureUnionize,rma::criteria,"
           f"[section_data_set_id$eq{dataset_id}],"
           f"[structure_id$in{id_str}]"
           f"&num_rows=50")
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    rows = r.json()["msg"]
    # map structure_id → expression_energy
    id2acr = {v: k for k, v in structure_ids.items()}
    return {id2acr[row["structure_id"]]: row["expression_energy"]
            for row in rows if row["structure_id"] in id2acr}


# ── Main ──────────────────────────────────────────────────────────────────────
print("Fetching structure IDs...")
struct_ids = get_structure_ids(STRUCTURE_ACRONYMS)
print(f"  Found: {struct_ids}")

results = []   # list of {gene, structure, expression_energy}

for gene in GENES:
    print(f"  {gene}...", end=" ")
    datasets = get_ish_datasets(gene)
    if not datasets:
        print("no datasets found, skipping")
        continue
    # Use first (most representative) dataset
    expr = get_structure_expression(datasets[0], struct_ids)
    for acr in STRUCTURE_ACRONYMS:
        results.append({
            "gene": gene,
            "structure": acr,
            "structure_label": STRUCTURE_LABELS[acr].replace("\n", " "),
            "expression_energy": expr.get(acr, None),
        })
    print(f"dataset {datasets[0]}, {expr}")
    time.sleep(0.3)  # be polite to the API

# Write TSV
with open(OUT / "atlas_expression.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["gene", "structure", "structure_label",
                                       "expression_energy"], delimiter="\t")
    w.writeheader()
    w.writerows(results)

print(f"\nWrote {len(results)} rows to atlas_expression.tsv")
