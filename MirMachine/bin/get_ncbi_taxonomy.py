#!/usr/bin/env python

import ete3
from ete3 import NCBITaxa


def get_ncbi(update: bool = False) -> ete3.ncbi_taxonomy.ncbiquery.NCBITaxa:
    """
    Get an instance of ete3.ncbi_taxonomy.ncbiquery.NCBITaxa

    Parameters
    ----------
    update : bool
        Update the NCBI taxonomy database

    Returns
    -------
    ete3.ncbi_taxonomy.ncbiquery.NCBITaxa
        Instance of ete3.ncbi_taxonomy.ncbiquery.NCBITaxa
    """
    ncbi = NCBITaxa(dbfile="taxa.sqlite")

    if update:
        ncbi.update_taxonomy_database()
    return ncbi


if __name__ == "__main__":
    ncbi = get_ncbi(update=True)