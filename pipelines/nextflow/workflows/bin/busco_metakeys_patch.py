#!/usr/bin/env python3
import argparse
import json
import re
from pathlib import Path
from typing import Dict, Union

def parse_busco_file(file_path: str, db: str) -> Dict[str, Union[str, float]]:
    """
    Parses a BUSCO result file and extracts relevant data into a dictionary.

    Args:
        file_path (str): The path to the BUSCO result file.
        db(str): Core db name.

    Returns:
        Dict[str, str]: A dictionary containing parsed BUSCO data, including the dataset,
                        completeness values, and mode (proteins or genome).
    """
    

    # Declare the dictionary to accept str as keys and str or float as values
    data: Dict[str, Union[str, float]] = {}
    data["core_db"] = db
    # Open and read the file
    with open(file_path, "r") as file:
        content = file.read()

    # Extract the BUSCO version
    version_match = re.search(r"BUSCO version is: ((\d+\.\d+.\d+))", content)
    if version_match:
        data["genebuild.busco_version"] = version_match.group(1)

    # Extract the BUSCO dataset
    dataset_match = re.search(r"The lineage dataset is: ([\w_]+)", content)
    if dataset_match:
        data["genebuild.busco_dataset"] = dataset_match.group(1)

    # Extract the BUSCO mode
    mode_match = re.search(r"BUSCO was run in mode: ([\w_]+)", content).group(1)
    if mode_match in ("genome", "euk_genome_met", "euk_genome_min"):
        data["genebuild.busco_mode"] = "genome"
    elif mode_match == "proteins":
        data["genebuild.busco_mode"] = "protein"
    else:
        mode_match = None    

    # Extract the BUSCO summary line with completeness values
    if mode_match == "euk_genome_min":
        result_match = re.search(
            r"C:(\d+\.\d+)%\[S:(\d+\.\d+)%.*,D:(\d+\.\d+)%\],F:(\d+\.\d+)%.*,M:(\d+\.\d+)%,n:(\d+),E:(\d+\.\d+)%",
            content,
        )
    else:
        result_match = re.search(
            r"C:(\d+\.\d+)%\[S:(\d+\.\d+)%.*,D:(\d+\.\d+)%\],F:(\d+\.\d+)%.*,M:(\d+\.\d+)%,n:(\d+)", content
        )
    if result_match:
        completeness = result_match.group(1)
        single_copy = result_match.group(2)
        duplicated = result_match.group(3)
        fragmented = result_match.group(4)
        missing = result_match.group(5)
        total_buscos = result_match.group(6)
        if mode_match == "euk_genome_min":
            erroneus = result_match.group(7)
            # Store the BUSCO completeness summary with erroneous
            data[
                f'genebuild.busco_{data["genebuild.busco_mode"]}'
            ] = f"C:{completeness}%[S:{single_copy}%,D:{duplicated}%],F:{fragmented}%,M:{missing}%,n:{total_buscos},E:{erroneus}%"
        else:
            # Store the BUSCO completeness summary
            data[
                f'genebuild.busco_{data["genebuild.busco_mode"]}'
            ] = f"C:{completeness}%[S:{single_copy}%,D:{duplicated}%],F:{fragmented}%,M:{missing}%,n:{total_buscos}"

        # Unpack the BUSCO values into individual fields
        data[f'genebuild.busco_{data["genebuild.busco_mode"]}_completeness'] = float(completeness)
        data[f'genebuild.busco_{data["genebuild.busco_mode"]}_single_copy'] = float(single_copy)
        data[f'genebuild.busco_{data["genebuild.busco_mode"]}_duplicated'] = float(duplicated)
        data[f'genebuild.busco_{data["genebuild.busco_mode"]}_fragmented'] = float(fragmented)
        data[f'genebuild.busco_{data["genebuild.busco_mode"]}_missing'] = float(missing)
        data[f'genebuild.busco_{data["genebuild.busco_mode"]}_total'] = int(total_buscos)
        if mode_match == "euk_genome_min":
            data[f'genebuild.busco_{data["genebuild.busco_mode"]}_erroneus'] = float(erroneus)
    return data


# Function to generate SQL patches
def generate_sql_patches(db_name: str, json_data: Dict[str, Union[str, float]], species_id: int = 1, table_name: str = "meta"):
    sql_statements = []
    sql_statements.append(f"USE {db_name};\n")  # Replace with your actual DB name

    # Iterate through the JSON key-value pairs
    for key, value in json_data.items():
        if value is None:
            # Skip if the value is None (or can handle it differently if needed)
            continue
        # Convert value to string and escape single quotes if necessary
        value_str = str(value).replace("'", "''")
        # Create the SQL INSERT statement
        sql_statements.append(
            f"INSERT INTO {table_name} (species_id, meta_key, meta_value) VALUES ({species_id}, '{key}', '{value_str}');\n"
        )

    return "".join(sql_statements)


def main():
    """
    Main function to parse a BUSCO result file and output the parsed data in JSON format.

    It expects the file path to the BUSCO result as a command-line argument.
    """

    # Set up argument parser
    parser = argparse.ArgumentParser(description="Parse a BUSCO result file and generate JSON output.")
    parser.add_argument("-file", type=str, help="Path to the BUSCO result file")
    parser.add_argument("-db", type=str, help="Core db")
    parser.add_argument("-input_dir", type=str, help="Path for directory containing the busco output files")
    parser.add_argument("-output_dir", type=str, help="Path for output directory")
    # Parse arguments
    args = parser.parse_args()
    if args.file:
        # Parse the BUSCO file and generate the JSON
        busco_data = parse_busco_file(args.file, args.db)
        # Determine the file name based on the mode (protein or genome)
        busco_mode = busco_data.get("genebuild.busco_mode", "unknown")
        output_file_name = f"{args.db}_busco_{busco_mode}_metakey.json"

        # Convert the dictionary to a JSON object
        busco_json = json.dumps(busco_data, indent=4)

        # Write the JSON output to the dynamically nasmed file
        with open(Path(args.output_dir) / output_file_name, "w") as outfile:
            outfile.write(busco_json)

        # Output the JSON
        print(busco_json)
    elif args.input_dir:
        # Find all files that end with 'busco_short_summary'
        busco_files = list(Path(args.input_dir).rglob("*busco_short_summary.txt"))
        with open(Path(args.output_dir) / f"{args.db}.sql", "a") as f:
            for file in busco_files:
                print(file)
                busco_data = parse_busco_file(file, args.db)
                # Determine the file name based on the mode (protein or genome)
                busco_mode = busco_data.get("genebuild.busco_mode", "unknown")
                output_file_name = f"{args.db}_busco_{busco_mode}_metakey.json"

                # Convert the dictionary to a JSON object
                busco_json = json.dumps(busco_data, indent=4)

                # Write the JSON output to the dynamically named file
                with open(output_file_name, "w") as outfile:
                    outfile.write(busco_json)

                # Output the JSON
                print(busco_json)

                # Generate SQL patches from the JSON
                sql_patches = generate_sql_patches(args.db, busco_data)

                f.write(sql_patches)


if __name__ == "__main__":
    main()
