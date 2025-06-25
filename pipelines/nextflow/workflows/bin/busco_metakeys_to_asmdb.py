#!/usr/bin/env python3
# pylint: disable=missing-module-docstring
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import argparse
import json
import pymysql
import re
from pathlib import Path
from typing import Dict, Optional, Union
import os

def load_json(filepath):
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"{filepath} does not exist")
    with open(filepath, 'r') as f:
        return json.load(f)

def get_asm_metadata(gca: str, metadata_params) -> int:
    """
     Get assembly_id, required to insert data to assembly_metrics table
    """
    try:
        connection = pymysql.connect(**metadata_params)
        with connection:
            with connection.cursor() as cursor:
                assembly_id_query = f""" SELECT assembly_id FROM assembly
                WHERE CONCAT(gca_chain, '.', gca_version) = '{gca}';
                """
                cursor.execute(assembly_id_query)
                output = cursor.fetchone()
                if not output:
                    raise ValueError(f"No assembly_id found for GCA: {gca}")
                print(f"Assembly id for {gca} is {output[0]}")
        return output[0]
    except Exception as e:
        print(f"Failed to fetch assembly id for {gca}: {e}")
        raise

def parse_busco_file(file_path: str) -> Tuple[Dict[str, Union[str, int]], str]:
    """
    Parses a BUSCO result file and extracts relevant data into a dictionary.

    Args:
        file_path (str): The path to the BUSCO result file.

    Returns:
        Dict[str, str]: A dictionary containing parsed BUSCO data, including the dataset,
                        completeness values, and mode (proteins or genome).
        str: a string specifiying busco mode
    """

    # Declare the dictionary to accept str as keys and str or float as values
    data: Dict[str, Union[str, int]] = {}
    # Open and read the file
    with open(file_path, "r") as file:
        content = file.read()

    # Define regular expressions to match the relevant numbers
    version_pattern : Optional[re.Match[str]] = re.search(r"BUSCO version is: ((\d+\.\d+.\d+))", content)
    dataset_pattern : Optional[re.Match[str]] = re.search(r"The lineage dataset is: ([\w_]+)", content)
    mode_pattern: Optional[re.Match[str]] = re.search(r"BUSCO was run in mode: ([\w_]+)", content)
    completeness_pattern: Optional[re.Match[str]] = re.search(r"(\d+)\s+Complete BUSCOs \(C\)", content)
    single_copy_pattern: Optional[re.Match[str]] = re.search(
        r"(\d+)\s+Complete and single-copy BUSCOs \(S\)", content
    )
    duplicates_pattern: Optional[re.Match[str]] = re.search(
        r"(\d+)\s+Complete and duplicated BUSCOs \(D\)", content
    )
    fragmented_pattern: Optional[re.Match[str]] = re.search(r"(\d+)\s+Fragmented BUSCOs \(F\)", content)
    missing_pattern: Optional[re.Match[str]] = re.search(r"(\d+)\s+Missing BUSCOs \(M\)", content)

    # Initialize mode_match as None or str
    mode_match: Optional[str] = None

    # If match is not None, extract the group and assign it to mode_match
    if mode_pattern is not None:
        mode_match = mode_pattern.group(1)
    if mode_match in ("genome", "euk_genome_met", "euk_genome_min"):
        busco_mode = "genome"
    elif mode_match == "proteins":
        busco_mode = "protein"
    else:
        mode_match = None

    version = str(version_pattern.group(1)) if version_pattern else None
    dataset = str(dataset_pattern.group(1)) if dataset_pattern else None
    completeness = int(completeness_pattern.group(1)) if completeness_pattern else None
    single_copy = int(single_copy_pattern.group(1)) if single_copy_pattern else None
    duplicated = int(duplicates_pattern.group(1)) if duplicates_pattern else None
    fragmented = int(fragmented_pattern.group(1)) if fragmented_pattern else None
    missing = int(missing_pattern.group(1)) if missing_pattern else None

    # Extract the BUSCO summary line with completeness values
    if mode_match == "euk_genome_min":
        score_match = re.search(
            r"C:(\d+\.\d+)%\[S:(\d+\.\d+)%.*,D:(\d+\.\d+)%\],F:(\d+\.\d+)%.*,M:(\d+\.\d+)%,n:(\d+),E:(\d+\.\d+)%",  # pylint: disable=line-too-long
            content,  # pylint: disable=line-too-long
        )
    else:
        score_match = re.search(
            r"C:(\d+\.\d+)%\[S:(\d+\.\d+)%.*,D:(\d+\.\d+)%\],F:(\d+\.\d+)%.*,M:(\d+\.\d+)%,n:(\d+)", content
        )

    if score_match:
        score = score_match.group(0)
        total_buscos = score_match.group(6)
        if mode_match == "euk_genome_min":
            erroneus = score_match.group(7)

        if mode_match in ("genome", "euk_genome_met", "euk_genome_min"):
            # Extract the BUSCO version
            data["assembly.busco_version"] = str(version)
            # Extract the BUSCO dataset
            data["assembly.busco_dataset"] = str(dataset)
            # Store the BUSCO completeness summary with erroneous
            data["assembly.busco"] = str(score)
            data["assembly.busco_mode"] = busco_mode
            # Store the BUSCO values into individual fields
            data["assembly.busco_completeness"] = str(completeness)
            data["assembly.busco_single_copy"] = str(single_copy)
            data["assembly.busco_duplicated"] = str(duplicated)
            data["assembly.busco_fragmented"] = str(fragmented)
            data["assembly.busco_missing"] = str(missing)
            data["assembly.busco_total"] = int(total_buscos)
            if mode_match == "euk_genome_min":
                data["assembly.busco_erroneus"] = str(erroneus)

            data["assembly.busco"] = str(score)  # pylint: disable=line-too-long

        else:
            # Extract the BUSCO version
            data["genebuild.busco_version"] = str(version)
            # Extract the BUSCO dataset
            data["genebuild.busco_dataset"] = str(dataset)
            # Store the BUSCO completeness summary
            data["genebuild.busco"] = str(score)
            data["genebuild.busco_mode"] = busco_mode
            # Store the BUSCO values into individual fields
            data["genebuild.busco_completeness"] = str(completeness)
            data["genebuild.busco_single_copy"] = str(single_copy)
            data["genebuild.busco_duplicated"] = str(duplicated)
            data["genebuild.busco_fragmented"] = str(fragmented)
            data["genebuild.busco_missing"] = str(missing)
            data["genebuild.busco_total"] = int(total_buscos)
    #print(f"Output from parse_busco_file. data: {data}")
    return data, busco_mode

def execute_query(busco_ncbi, assembly_id, metadata_params):
    """
    Inserts BUSCO NCBI into the assembly_metrics table for a given assembly_id.
    """
    try:
        connection = pymysql.connect(**metadata_params)
        with connection:
            with connection.cursor() as cursor:
                for key, value in busco_ncbi.items():
                    insert_query = """
                    INSERT INTO assembly_metrics (assembly_id, metrics_name, metrics_value)
                    VALUES (%s, %s, %s)
                    ON DUPLICATE KEY UPDATE metrics_value = VALUES(metrics_value);
                    """
                    cursor.execute(insert_query, (assembly_id, key, value))
    except Exception as e:
        print(f"Failed to insert query {(insert_query, (assembly_id, key, value))}: {e}")
        raise

def main():
    """
    Module to parse BUSCO genome (NCBI) output file to be store in the assembly metadata DB (assembly metrics table).

    It expects the file path to the BUSCO output file, full accession (GCA) and DB params.
    """

    # Set up argument parser
    parser = argparse.ArgumentParser(description="Parse a BUSCO result file and generate JSON output.")
    parser.add_argument("-file", type=str, help="Path to the BUSCO result file")
    parser.add_argument("-gca", type=str, help="GCA accession")
    parser.add_argument("-output_dir", type=str, help="Path for output directory")
    parser.add_argument("-asm_metadata", type=str, help="Path to json file with parms to assembly metadata db")

    # Parse arguments
    args = parser.parse_args()

    # Load assembly DB params
    metadata_params = load_json(args.asm_metadata)

    # Process the busco file and get a dictionary
    busco_ncbi,busco_mode = parse_busco_file(args.file)

    # Get assembly id, require to insert to assembly_metics
    assembly_id = get_asm_metadata(args.gca, metadata_params) 
    
    # Execute query 
    execute_query(busco_ncbi, assembly_id, metadata_params)

if __name__ == "__main__":
    main()
