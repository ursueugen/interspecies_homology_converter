import shutil
from pathlib import Path
import pandas as pd
from pybiomart import Dataset


FILEPATH = "genage_Celegans_allgenes_proantilongevity.tsv"
RESULTS_DIR = Path("Celegans_to_Dmelanogaster_Homology")
RESULTPATH = RESULTS_DIR / (FILEPATH.split(".")[0] 
                            + "_with_Homologs.tsv")


ENTREZ_ID_COL = "Entrez Gene ID"

LOOKUP_FILENAME = RESULTS_DIR / "CELEGANS_DMELANOGASTER_HOMOLOGY_LOOKUP.tsv"

HOST = "http://www.ensembl.org"
CELEGANS_DATASET_NAME = "celegans_gene_ensembl"
DROSO_DATASET_NAME = "dmelanogaster_gene_ensembl"

ENSEMBL_ID_ATTRIBUTE = "ensembl_gene_id"
ENTREZ_ID_ATTRIBUTE = "entrezgene_id"

DROSO_HOMO_GENE = "dmelanogaster_homolog_ensembl_gene"
DROSO_HOMO_GENE_NAME = "dmelanogaster_homolog_associated_gene_name"
DROSO_HOMO_HOMO_TYPE = "dmelanogaster_homolog_orthology_type"
DROSO_HOMO_CONFIDENCE = "dmelanogaster_homolog_orthology_confidence"
DROSO_HOMO_ATTRIBUTES = [DROSO_HOMO_GENE, 
                         DROSO_HOMO_GENE_NAME,
                         DROSO_HOMO_HOMO_TYPE,
                         DROSO_HOMO_CONFIDENCE]


def main():
    
    if RESULTS_DIR.is_dir():
        shutil.rmtree(RESULTS_DIR)
    
    RESULTS_DIR.mkdir()
    
    lookup = get_homology_lookup()
    lookup = add_entrez_ids(lookup)
    
    
    df = read(FILEPATH)
#     df.dropna(subset=[ENTREZ_ID_COL], axis=0, inplace=True)
    
    df_final = add_orthology_to_dataframe(df, lookup)
    
    df_final.to_csv(RESULTPATH, header=True, index=False)
    

def read(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t')
    return df


def get_homology_lookup() -> pd.DataFrame:
    """
    Returns lookup table consisting of ensembl id of reference species (C.elegans)
        and ensembl id, gene symbol, orthology type and orthology confidence
        of the other species (D.melanogaster).
    """
    
    dataset = Dataset(name=CELEGANS_DATASET_NAME,
                     host=HOST)
    
    attributes = [ENSEMBL_ID_ATTRIBUTE] + DROSO_HOMO_ATTRIBUTES
    df_lookup = dataset.query(attributes=attributes,
                    filters=None)

    df_lookup.to_csv(LOOKUP_FILENAME, header=True, index=True)
    
    return df_lookup


def get_species_ens_entrez_lookup(dataset_name: str) -> pd.DataFrame:
    """
    Returns lookup table for a ensembl dataset name with 2 columns:
        ensembl id, entrez id.
    """
    
    dataset = Dataset(name=dataset_name,
                     host=HOST)
    
    df_lookup = dataset.query(attributes=[
        ENSEMBL_ID_ATTRIBUTE,
        ENTREZ_ID_ATTRIBUTE],
        filters=None)
    
    df_lookup.to_csv(RESULTS_DIR / ("{}_ENS_ENTREZ_LOOKUP_.csv"
                     .format(dataset_name.split("_")[0].upper())),
                     header=True, index=True)
    
    return df_lookup


def add_entrez_ids(lookup: pd.DataFrame) -> pd.DataFrame:
    """
    Adds entrez ids for celeg and dmelanogaster to a lookup table.
    """
    
    celeg_ens2entrez = get_species_ens_entrez_lookup(CELEGANS_DATASET_NAME)
    celeg_ens2entrez.columns = ['celeg_ensembl_id',
                                 'celeg_entrez_id']
    
    droso_ens2entrez = get_species_ens_entrez_lookup(DROSO_DATASET_NAME)
    droso_ens2entrez.columns = ['dmelanogaster_ensembl_id',
                                'dmelanogaster_entrez_id']
    
    lookup_with_entrez = pd.merge(lookup, celeg_ens2entrez,
                                 left_on="Gene stable ID",
                                 right_on="celeg_ensembl_id",
                                 how="left")
    
    lookup_with_entrez = pd.merge(lookup_with_entrez, droso_ens2entrez,
                                 left_on="Drosophila melanogaster gene stable ID",
                                 right_on="dmelanogaster_ensembl_id",
                                 how="left")
    
    lookup_with_entrez.to_csv(LOOKUP_FILENAME, header=True, index=False)
    return lookup_with_entrez


def add_orthology_to_dataframe(df, lookup):
    """
    Adds orthology data to input dataframe containing a column with entrez ids
        of C.elegans genes.
    """
    
    mask_na = df[ENTREZ_ID_COL].isna()
    
    df_final = pd.merge(df.loc[~mask_na], lookup,
                       left_on=ENTREZ_ID_COL,
                       right_on="celeg_entrez_id",
                       how="left")
    
    df_final = pd.concat([df_final,
                        df.loc[mask_na]],
                        axis=0,
                        sort=False)
    
    return df_final


if __name__ == "__main__":
    main()