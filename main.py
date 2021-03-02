# -*- coding: utf-8 -*-

"""Build the "Not Scared of Chemistry" Knowledge Graph."""

import lzma
import zipfile
from typing import TextIO

import click
import more_click
from tqdm import tqdm

import bioversions
import pystow

# unlike BioGRID, ExCAPE-DB can be considered a static resource
EXCAPE_URL = 'https://zenodo.org/record/2543724/files/pubchem.chembl.dataset4publication_inchi_smiles_v2.tsv.xz'
EXCAPE_VERSION = 'v2'


def _excape(file: TextIO, human_only: bool = False) -> None:
    """Pre-process ExCAPE-DB.

    ExCAPE-DB is a database of chemical modulations of proteins built as a curated subset of
    ChEBML and PubChem

    Future directions:
    - Add a variable pXC50 cutoff besides 6.0
    """
    excape_path = pystow.ensure('bio2bel', 'excapedb', EXCAPE_VERSION, url=EXCAPE_URL)
    with lzma.open(excape_path, mode='rt') as infile:
        _header = next(infile)
        for i, line in enumerate(tqdm(infile, unit_scale=True, desc='ExCAPE-DB')):
            line = line.strip().split('\t')
            if human_only and line[7] != '9606':  # Taxonomy ID must be human
                continue
            if line[3] != 'A':  # Activity_Flag is active (A) instead of not active (N)
                continue
            target = line[2]
            try:
                int(target)
            except ValueError:
                tqdm.write(f'failure on line {i}')
                continue
            else:
                print(f'inchikey:{line[0]}', 'modulates', f'ncbigene:{line[2]}', sep='\t', file=file)


def _biogrid(file: TextIO, version: str, human_only: bool = False) -> None:
    """Pre-process the given version of BioGRID.

    BioGRID is a manually curated database of protein-protein and protein-complex interactions.
    """
    url = f'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/' \
          f'BIOGRID-{version}/BIOGRID-ALL-{version}.tab3.zip'
    path = pystow.ensure('bio2bel', 'biogrid', version, url=url)

    with zipfile.ZipFile(path) as zip_file:
        with zip_file.open(f'BIOGRID-ALL-{version}.tab3.txt') as infile:
            lines = (
                line.decode('utf-8').strip().split('\t')
                for line in tqdm(infile, unit_scale=True, desc=f'BioGRID v{version}')
            )

            header = next(lines)
            header_dict = {entry: i for i, entry in enumerate(header)}

            source_key = header_dict['Entrez Gene Interactor A']
            target_key = header_dict['Entrez Gene Interactor B']
            organism_a_key = header_dict['Organism Name Interactor A']
            organism_b_key = header_dict['Organism Name Interactor B']

            for line in lines:
                if human_only and (line[organism_a_key] != 'Homo sapiens' or line[organism_b_key] != 'Homo sapiens'):
                    continue
                print(
                    f'ncbigene:{line[source_key]}',
                    'interacts',
                    f'ncbigene:{line[target_key]}',
                    sep='\t',
                    file=file,
                )


def _homologene(file: TextIO, version: str) -> None:
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
    url = f'https://ftp.ncbi.nih.gov/pub/HomoloGene/build{version}/homologene.data'
    df = pystow.ensure_csv('bio2bel', 'homologene', version, url=url, read_csv_kwargs={
        'usecols': [0, 2],
        'sep': '\t',
    })
    for homologene_id, ncbigene_id in tqdm(df.values, unit_scale=True, desc=f'HomoloGene v{version}'):
        print(ncbigene_id, 'homologyGroup', homologene_id, sep='\t', file=file)


@click.command()
@more_click.verbose_option
def main():
    """Build NSoC-KG."""
    excape_cut_path = pystow.get('nsockg', 'excapedb', EXCAPE_VERSION, 'excape.tsv')
    with open(excape_cut_path, 'w') as file:
        _excape(file)

    biogrid_version = bioversions.get_version('biogrid')
    biogrid_cut_path = pystow.get('nsockg', 'biogrid', biogrid_version, 'biogrid_cut.tsv')
    with open(biogrid_cut_path, 'w') as file:
        _biogrid(file, biogrid_version)

    homolgene_version = bioversions.get_version('homologene')
    homolgene_cut_path = pystow.get('nsockg', 'homologene', homolgene_version, 'homologene.tsv')
    with open(homolgene_cut_path, 'w') as file:
        _homologene(file, homolgene_version)


if __name__ == '__main__':
    main()
