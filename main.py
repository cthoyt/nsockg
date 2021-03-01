# -*- coding: utf-8 -*-

"""Build the "Not Scared of Chemistry" Knowledge Graph."""

import lzma
import zipfile

import click
import more_click
from tqdm import tqdm

import bioversions
import pystow

# unlike BioGRID, ExCAPE-DB can be considered a static resource
EXCAPE_URL = 'https://zenodo.org/record/2543724/files/pubchem.chembl.dataset4publication_inchi_smiles_v2.tsv.xz'
EXCAPE_VERSION = 'v2'


def _excape(file) -> None:
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
            if line[7] != '9606':  # Taxonomy ID must be human
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


def _biogrid(file, version):
    """BioGRID

    BioGRID is a manually curated database of protein-protein and protein-complex interactions

    :param file:
    :return:
    """
    biogrid_url = f'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/' \
                  f'BIOGRID-{version}/BIOGRID-ALL-{version}.tab3.zip'
    biogrid_path = pystow.ensure('bio2bel', 'biogrid', version, url=biogrid_url)

    with zipfile.ZipFile(biogrid_path) as zip_file:
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
                if line[organism_a_key] != 'Homo sapiens':
                    continue
                if line[organism_b_key] != 'Homo sapiens':
                    continue
                print(
                    f'ncbigene:{line[source_key]}',
                    'interacts',
                    f'ncbigene:{line[target_key]}',
                    sep='\t',
                    file=file,
                )


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


if __name__ == '__main__':
    main()
