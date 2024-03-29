# -*- coding: utf-8 -*-

"""Build the "Not Scared of Chemistry" Knowledge Graph.

.. seealso:: https://doi.org/10.5281/zenodo.4574554
"""

import datetime
import getpass
import json
from typing import Iterable, Sequence, TextIO

import bioversions
import click
import more_click
import pystow
from tabulate import tabulate
from tqdm import tqdm
from zenodo_client import Creator, Metadata, ensure_zenodo

# unlike BioGRID, ExCAPE-DB can be considered a static resource
EXCAPE_URL = "https://zenodo.org/record/2543724/files/pubchem.chembl.dataset4publication_inchi_smiles_v2.tsv.xz"
EXCAPE_VERSION = "v2"

DISGENET_URL = "https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz"

NSOCKG_MODULE = pystow.module("nsockg")
BIO = pystow.module("bio")

metadata = Metadata(
    title="Not Scared of Chemistry Knowledge Graph",
    upload_type="dataset",
    description="A combination of ExCAPE-DB, BioGRID, HomoloGene, and chemical similarities in a knowledge graph.",
    creators=[
        Creator(
            name="Hoyt, Charles Tapley",
            affiliation="Harvard Medical School",
            orcid="0000-0003-4423-4370",
        ),
    ],
)


# EXCAPE - CC-BY-SA-4.0 License
# BioGRID - MIT License
# HomoloGene - ??
# DisGeNet - Attribution-NonCommercial-ShareAlike 4.0 International License


@click.command()
@more_click.verbose_option
def main():
    """Build NSoC-KG."""
    biogrid_version = bioversions.get_version("biogrid")
    homolgene_version = bioversions.get_version("homologene")
    disgenet_version = bioversions.get_version("disgenet")
    excape_version = EXCAPE_VERSION
    versions = {
        "biogrid": biogrid_version,
        "homologene": homolgene_version,
        "excape": excape_version,
        "disgenet": disgenet_version,
    }

    statistics = {}
    triples_path = NSOCKG_MODULE.join(name="triples.tsv")
    with triples_path.open("w") as file:
        _excape(statistics, file, excape_version)
        _biogrid(statistics, file, biogrid_version)
        _homologene(statistics, file, homolgene_version)
        _disgenet(statistics, file, disgenet_version)

    # Count everything
    statistics["total"] = sum(statistics.values())

    rows = [(key, versions[key], statistics[key]) for key in sorted(versions)]
    rows.append(("total", "", statistics["total"]))
    print(tabulate(rows, headers=["Source", "Version", "Edges"]))

    metadata_path = NSOCKG_MODULE.join(name="metadata.json")
    with metadata_path.open("w") as file:
        json.dump(
            fp=file,
            indent=2,
            obj={
                "date": datetime.datetime.now().strftime("%Y-%m-%d"),
                "exporter": getpass.getuser(),
                "versions": versions,
                "statistics": statistics,
            },
        )

    # Automatically upload this revision to Zenodo
    ensure_zenodo(
        key="nsockg",
        data=metadata,
        paths=[
            triples_path,
            metadata_path,
        ],
    )


def _disgenet(statistics, file: TextIO, version: str) -> None:
    module = BIO.submodule("disgenet", version)
    df = module.ensure_csv(
        url=DISGENET_URL,
        read_csv_kwargs={"dtype": {"geneId": str}},
    )
    count = 0
    for index, row in tqdm(
        df.iterrows(), total=len(df.index), unit_scale=True, desc=f"DisGeNet v{version}"
    ):
        ncbigene_id = row["geneId"].strip()
        disease_umls_id = row["diseaseId"]
        count += 1
        print(
            f"ncbigene:{ncbigene_id}",
            "associated",
            f"umls:{disease_umls_id}",
            sep="\t",
            file=file,
        )
    statistics["disgenet"] = count


def _excape(
    statistics,
    file: TextIO,
    version: str,
    human_only: bool = False,
    url: str = EXCAPE_URL,
) -> None:
    """Pre-process ExCAPE-DB.

    ExCAPE-DB is a database of chemical modulations of proteins built as a curated subset of
    ChEBML and PubChem

    Future directions:
    - Add a variable pXC50 cutoff besides 6.0
    """
    module = BIO.submodule("excapedb", version)
    with module.ensure_open_lzma(url=url) as infile:
        _header = next(infile)
        it = tqdm(infile, unit_scale=True, desc=f"ExCAPE-DB {version}")
        for i, line in enumerate(it):
            line = line.strip().split("\t")
            if human_only and line[7] != "9606":  # Taxonomy ID must be human
                continue
            if line[3] != "A":  # Activity_Flag is active (A) instead of not active (N)
                continue
            target = line[2]
            try:
                int(target)
            except ValueError:
                it.write(f"failure on line {i}")
                continue
            else:
                print(
                    f"inchikey:{line[0]}",
                    "modulates",
                    f"ncbigene:{line[2]}",
                    sep="\t",
                    file=file,
                )

        statistics["excape"] = i


def _biogrid(statistics, file: TextIO, version: str, human_only: bool = False) -> None:
    """Pre-process the given version of BioGRID.

    BioGRID is a manually curated database of protein-protein and protein-complex interactions.
    """
    url = (
        f"https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/"
        f"BIOGRID-{version}/BIOGRID-ALL-{version}.tab3.zip"
    )
    inner_path = f"BIOGRID-ALL-{version}.tab3.txt"
    module = BIO.submodule("biogrid", version)
    with module.ensure_open_zip(url=url, inner_path=inner_path) as infile:
        lines = (
            line.decode("utf-8").strip().split("\t")
            for line in tqdm(infile, unit_scale=True, desc=f"BioGRID v{version}")
        )

        header = next(lines)
        header_dict = {entry: i for i, entry in enumerate(header)}

        source_key = header_dict["Entrez Gene Interactor A"]
        target_key = header_dict["Entrez Gene Interactor B"]
        organism_a_key = header_dict["Organism Name Interactor A"]
        organism_b_key = header_dict["Organism Name Interactor B"]

        count = 0
        for line in lines:
            if human_only and (
                line[organism_a_key] != "Homo sapiens"
                or line[organism_b_key] != "Homo sapiens"
            ):
                continue

            count += 1
            print(
                f"ncbigene:{line[source_key]}",
                "interacts",
                f"ncbigene:{line[target_key]}",
                sep="\t",
                file=file,
            )
        statistics["biogrid"] = count


def _homologene(statistics, file: TextIO, version: str) -> None:
    """Pre-process the orthology data from HomoloGene.

    :param file:
    :param version:
    :return:

    The README at https://ftp.ncbi.nih.gov/pub/HomoloGene/README states
    that HomoloGene has the following data:

    1) HID (HomoloGene group id)
    2) Taxonomy ID
    3) Gene ID
    4) Gene Symbol
    5) Protein gi
    6) Protein accession

    """
    url = f"https://ftp.ncbi.nih.gov/pub/HomoloGene/build{version}/homologene.data"
    module = BIO.submodule("homologene", version)
    df = module.ensure_csv(
        url=url,
        read_csv_kwargs={
            "usecols": [0, 2],
        },
    )
    count = 0
    it = tqdm(df.values, unit_scale=True, desc=f"HomoloGene v{version}")
    for homologene_id, ncbigene_id in it:
        count += 1
        print(
            f"ncbigene:{ncbigene_id}",
            "homologyGroup",
            f"homologene:{homologene_id}",
            sep="\t",
            file=file,
        )
    statistics["homologene"] = count


def cut(
    file: Iterable[str], sep: str, columns: Sequence[int]
) -> Iterable[Sequence[str]]:
    for line in file:
        line = line.strip().split(sep)
        yield tuple(line[column] for column in columns)


if __name__ == "__main__":
    main()
