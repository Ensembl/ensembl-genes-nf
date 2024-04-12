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
from pathlib import Path
from typing import Any, List
import requests


def get_dataset_match(ncbi_url: str, dataset: list) -> List[Any]:
    """
    Get taxonomy tree from ncbi taxonomy datasets and find the closest match with the input list


    Args:
        ncbi_url (str): Ncbi dataset url
        dataset (list): list of data to match

    Returns:
        str: closest match in the dataset list

    Raises:
        requests.HTTPError: If an HTTP error occurs during the API request.
        Exception: If any other error occurs during the function's operation.

    """

    try:
        response = requests.get(ncbi_url, timeout=10)
        response.raise_for_status()  # Raises an HTTPError if the response status code is 4XX or 5XX
        data = response.json()
        # data = json.loads(json_response)

        # Extract classification names
        classification = data["reports"][0]["taxonomy"]["classification"]
        classification_names = [value["name"] for key, value in classification.items()]
        match = []
        for i in reversed(classification_names):
            for l in dataset:
                if l.strip().lower() == i.strip().lower():
                    match.append(l.strip())
                    break  # Break out of the loop once an exact match is found
                if l.strip().startswith(str(i[: len(i) - 2].lower())):
                    match.append(l)
                    break  # Break out of the loop once a partial match is found
    except requests.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err}")
    except requests.exceptions.RequestException as req_err:
        print(f"Request error occurred: {req_err}")

    return match


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Clade selector arguments")
    parser.add_argument(
        "-d",
        "--datasets",
        type=str,
        help="Path to file containing list of datasets (one per line)",
        required=True,
    )
    parser.add_argument("-t", "--taxon_id", type=str, help="Taxon id ", required=True)
    parser.add_argument("--output", type=str, help="Output file", default="stdout")
    parser.add_argument(
        "--ncbi_url",
        type=str,
        help="NCBI dataset url",
        default="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/",
    )
    return parser.parse_args()


def main():
    """Entry-point."""
    args = parse_args()

    ncbi_url = f"{args.ncbi_url}/{args.taxon_id}/dataset_report"

    with open(Path(args.datasets), "r") as file:
        datasets = [line[: max(line.find(" "), 0) or None] for line in file]

    clade_match = get_dataset_match(ncbi_url, datasets)

    if not clade_match:
        raise ValueError("No match found")

    if args.output == "stdout":  # pylint:disable=no-else-return
        return clade_match[0].strip("\n")
    else:
        with open(args.output, "w+") as output:
            if clade_match[0] == args.species:
                output.write(clade_match[1])
            else:
                output.write(clade_match[0])

    return None


if __name__ == "__main__":
    main()
