"""
Script that given a list of clades (e.g. from the NCBI taxonomy) and a species name (e.g. from the Ensembl API) will return the closest match

Usage:
    match_clade.py --clades <clades> --species <species> [--output <output>]

Options:
    --clades <clades>       Path to file containing list of clades (one per line)
    --species <species>     Species name    
    --output <output>       Output file [default: ./clade_match.txt]

Example:
    match_clade.py --clades clades.txt --species "Homo sapiens" --output clade_match.txt
"""

import argparse
import ete3
from ete3 import NCBITaxa

try:
    import taxoniq
except ImportError:
    pass


def check_lineage_order(species: str, lineage: list) -> bool:
    """
    Return True if the lineage is in the correct order (from higher to lower)

    Parameters
    ----------
    lineage : list
        List of lineages

    Returns
    -------
    bool
        True if the lineage is in the correct order (from higher to lower)
    """
    if lineage[-1] == species:
        return True
    else:
        return False


def taxoniq_match(species: str, taxid: int, tool_lineage: list) -> str:
    """
    Check for a suitable match using taxoniq

    Parameters
    ----------
    species : str
        Species name
    taxid : int
        Taxonomy ID
    tool_lineage : list
        List of lineages

    Returns
    -------
    str
        Closest match
    """
    try:
        t = taxoniq.Taxon(taxid)
    except Exception as e:
        return e
    lineage = [t.scientific_name for t in t.ranked_lineage]
    if not check_lineage_order(t.scientific_name, lineage):
        lineage = lineage[::-1]

    if tool_lineage is None:
        raise ValueError("Tool lineages must be provided")

    matches = []
    for i in reversed(lineage):
        for l in tool_lineage:
            if l.strip().lower() == i.strip().lower():
                matches.append(l.strip())
            elif l.strip().startswith(str(i.strip()[: len(i) - 2].lower())):
                matches.append(l)
    return matches


def custom_get_name_translator(ncbi, names):
    """
    Custom implementation of get_name_translator to avoid using the problematic ete3 function.
    Given a list of taxid scientific names, returns a dictionary translating them into their corresponding taxids.
    
    Parameters
    ----------
    ncbi : ete3.ncbi_taxonomy.ncbiquery.NCBITaxa
        Instance of ete3.ncbi_taxonomy.ncbiquery.NCBITaxa
    names : list
        List of species names to translate
        
    Returns
    -------
    dict
        Dictionary mapping species names to taxids
    """
    name2id = {}
    name2origname = {}
    
    # Normalize names for case-insensitive comparison
    for n in names:
        name2origname[n.lower()] = n
    
    names = set(name2origname.keys())
    
    # SQL query needs proper escaping of single quotes
    query_names = []
    for name in name2origname.keys():
        # Replace single quotes with double single quotes for SQL
        escaped_name = name.replace("'", "''")
        query_names.append(f"'{escaped_name}'")
    
    query = ','.join(query_names)
    
    # First try the species table
    cmd = f"select spname, taxid from species where spname IN ({query})"
    result = ncbi.db.execute(cmd)
    
    for sp, taxid in result.fetchall():
        oname = name2origname[sp.lower()]
        if oname not in name2id:
            name2id[oname] = []
        name2id[oname].append(taxid)
    
    # Check for any names that weren't found
    missing = names - set([n.lower() for n in name2id.keys()])
    
    if missing:
        # Try the synonym table for missing names
        query_names = []
        for name in missing:
            escaped_name = name.replace("'", "''")
            query_names.append(f"'{escaped_name}'")
        
        query = ','.join(query_names)
        cmd = f"select spname, taxid from synonym where spname IN ({query})"
        result = ncbi.db.execute(cmd)
        
        for sp, taxid in result.fetchall():
            oname = name2origname[sp.lower()]
            if oname not in name2id:
                name2id[oname] = []
            name2id[oname].append(taxid)
    
    return name2id


def custom_get_taxid_translator(ncbi, taxids, try_synonyms=True):
    """
    Custom implementation of get_taxid_translator to avoid using the problematic ete3 function.
    Given a list of taxids, returns a dictionary with their corresponding scientific names.
    
    Parameters
    ----------
    ncbi : ete3.ncbi_taxonomy.ncbiquery.NCBITaxa
        Instance of ete3.ncbi_taxonomy.ncbiquery.NCBITaxa
    taxids : list
        List of taxids to translate
    try_synonyms : bool, optional
        Whether to try looking up taxids in the merged table, by default True
        
    Returns
    -------
    dict
        Dictionary mapping taxids to scientific names
    """
    # Convert everything to integers and remove None/empty values
    all_ids = set(map(int, taxids))
    all_ids.discard(None)
    all_ids.discard("")
    
    # Create the query string - no need to escape integers
    query = ','.join([f"'{v}'" for v in all_ids])
    
    # Query the species table
    cmd = f"select taxid, spname FROM species WHERE taxid IN ({query});"
    result = ncbi.db.execute(cmd)
    
    id2name = {}
    for tax, spname in result.fetchall():
        id2name[int(tax)] = spname

    # If we didn't find all taxids and try_synonyms is True, check the merged table
    if len(all_ids) != len(id2name) and try_synonyms:
        # Find taxids that weren't translated
        not_found_taxids = all_ids - set(id2name.keys())
        
        # Custom implementation of _translate_merged
        conv_all_taxids = set(not_found_taxids)
        cmd = f"select taxid_old, taxid_new FROM merged WHERE taxid_old IN ({','.join([f'{v}' for v in not_found_taxids])})"
        result = ncbi.db.execute(cmd)
        
        # Create mappings from old to new taxids
        old2new = {}
        new2old = {}
        for old, new in result.fetchall():
            old = int(old)
            new = int(new)
            conv_all_taxids.discard(old)
            conv_all_taxids.add(new)
            old2new[old] = new
            new2old[new] = old

        # If we found any merged taxids, look up their new names
        if old2new:
            query = ','.join([f"'{v}'" for v in new2old])
            cmd = f"select taxid, spname FROM species WHERE taxid IN ({query});"
            result = ncbi.db.execute(cmd)
            
            # Map the new taxid names back to the original old taxids
            for tax, spname in result.fetchall():
                old_taxid = new2old[int(tax)]
                id2name[old_taxid] = spname

    return id2name


def ete3_match(
    ncbi: ete3.ncbi_taxonomy.ncbiquery.NCBITaxa, species: str, tool_lineage: list
) -> str:
    """
    Check for a suitable match using ete3

    Parameters
    ----------
    ncbi : ete3.ncbi_taxonomy.ncbiquery.NCBITaxa
        Instance of ete3.ncbi_taxonomy.ncbiquery.NCBITaxa
    species : str
        Species name
    tool_lineage : list
        List of lineages

    Returns
    -------
    str
        Closest match
    """

    taxon_id = get_taxid(ncbi, species)
    lineage = ncbi.get_lineage(taxon_id)
    # names = ncbi.get_taxid_translator(lineage)
    names = custom_get_taxid_translator(ncbi, lineage)

    lineage_list = [names[taxid] for taxid in lineage]
    if not check_lineage_order(species, lineage_list):
        lineage_list = lineage_list[::-1]

    match = []
    for i in reversed(lineage_list):
        for l in tool_lineage:
            if l.strip().lower() == i.strip().lower():
                match.append(l.strip())
            elif l.strip().startswith(str(i[: len(i) - 2].lower())):
                match.append(l)

    return match


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
    ncbi = NCBITaxa()
    if update:
        ncbi.update_taxonomy_database()
    return ncbi


def get_taxid(ncbi: ete3.ncbi_taxonomy.ncbiquery.NCBITaxa, species_name: str) -> int:
    """
    Given a species name, return the taxonomy ID

    Parameters
    ----------
    species_name : str
        Species name

    Returns
    -------
    int
        Taxonomy ID
    """
    try:
        name2taxid = custom_get_name_translator(ncbi, [species_name])
        taxon_id = str(list(name2taxid.values())[0]).strip("[]")
        return taxon_id
    except IndexError:
        ncbi = get_ncbi(update=True)
        name2taxid = custom_get_name_translator(ncbi, [species_name])
        taxon_id = str(list(name2taxid.values())[0]).strip("[]")
        return taxon_id
    except:
        raise ValueError(ncbi.get_name_translator(["homo sapiens"]))
        # raise ValueError(f"Species name not found")


def get_clade_match(
    ncbi: ete3.ncbi_taxonomy.ncbiquery.NCBITaxa, clades: list, species_name: str
) -> str:
    """
    Given a list of clades and a species name, return the closest match

    Parameters
    ----------
    clades : list
        List of clades
    species_name : str
        Species name

    Returns
    -------
    str
        Closest match
    """
    species_name = " ".join(species_name.split(" ")[:2])
    if not check_lineage_order(species_name, clades):
        clades = clades[::-1]
    return ete3_match(ncbi, species_name, clades)


def is_deuterostome(
    ncbi: ete3.ncbi_taxonomy.ncbiquery.NCBITaxa, species_name: str
) -> str:
    """
    Determine if a species is a deuterostome or protostome.

    Parameters
    ----------
    ncbi : ete3.ncbi_taxonomy.ncbiquery.NCBITaxa
        Instance of ete3.ncbi_taxonomy.ncbiquery.NCBITaxa
    species_name : str
        Species name

    Returns
    -------
    str
        'Deuterostomia' or 'Protostomia'
    """
    taxid = get_taxid(ncbi, species_name)
    lineage = ncbi.get_lineage(taxid)
    names = custom_get_taxid_translator(ncbi, lineage)
    lineage_names = [names[taxid] for taxid in lineage]

    if "Deuterostomia" in lineage_names:
        return "deutero"
    elif "Protostomia" in lineage_names:
        return "proto"
    else:
        return "combined" # this is the default option


def main(args: argparse.Namespace) -> None:
    """
    Main function

    Parameters
    ----------
    args : argparse.Namespace
        Arguments from argparse
    """
    with open(args.clades, "r") as file:
        clades = [line[: max(line.find(" "), 0) or None] for line in file]

    ncbi = get_ncbi()
    clade_match = get_clade_match(ncbi, clades, args.species)
    deutero_proto = is_deuterostome(ncbi, args.species)

    if clade_match == []:
        raise ValueError("No match found")

    if args.output == "stdout":
        if clade_match[0] == args.species:
            print(f"{clade_match[1].strip()};{deutero_proto}", end="")
        else:
            print(f"{clade_match[0].strip()};{deutero_proto}", end="")
    else:
        with open(args.output, "w+") as output:
            if clade_match[0] == args.species:
                output.write(f"{clade_match[1]};{deutero_proto}")
            else:
                output.write(f"{clade_match[0]};{deutero_proto}")

    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--clades",
        type=str,
        help="Path to file containing list of clades (one per line)",
        required=True,
    )
    parser.add_argument("-s", "--species", type=str, help="Species name", required=True)
    parser.add_argument(
        "--output", type=str, help="Output file", default="./clade_match.txt"
    )
    args = parser.parse_args()
    main(args)
