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
import sys
from pathlib import Path
from typing import Any, Dict, Optional

import urllib.request
import urllib.error

import sys
import time


def get_dataset_match(ncbi_url: str, dataset: list, max_retries: int = 5, base_delay: int = 5):
    """
    Get taxonomy tree from NCBI taxonomy datasets and find the closest match
    with the input list. Retries on URL errors.
    """

    for attempt in range(1, max_retries + 1):
        matched_value = None  # ensure this always exists for each attempt

        try:
            with urllib.request.urlopen(ncbi_url, timeout=10) as response:
                data = response.read().decode('utf-8')
                json_data = json.loads(data)

                parents = json_data["reports"][0]["taxonomy"]["parents"]

                for parent_id in reversed(parents):
                    parent_id_str = str(parent_id)
                    if parent_id_str in dataset:
                        matched_value = dataset[parent_id_str]
                        break

            # If we got here without an exception, return (even if None)
            return matched_value

        except urllib.error.URLError as url_err:
            print(
                f"URL error occurred (attempt {attempt}/{max_retries}): {url_err}",
                file=sys.stderr,
            )
        except json.JSONDecodeError as json_err:
            print(
                f"Error decoding JSON from NCBI (attempt {attempt}/{max_retries}): {json_err}",
                file=sys.stderr,
            )

        # If not the last attempt, sleep with exponential backoff
        if attempt < max_retries:
            delay = base_delay * attempt
            time.sleep(delay)

    # All retries failed
    return None



def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Clade selector arguments")
    parser.add_argument(
        "-d",
        "--datasets",
        type=str,
        help="Path to JSON file containing BUSCO lineage datasets",
        required=True,
    )
    parser.add_argument(
        "-t",
        "--taxon_id",
        type=str,
        help="Taxon id",
        required=True,
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output file (default: stdout)",
        default="stdout",
    )
    parser.add_argument(
        "--ncbi_url",
        type=str,
        help="NCBI dataset base URL",
        default="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon",
    )
    return parser.parse_args()


def main():
    """Entry-point."""
    args = parse_args()

    # Load datasets JSON
    with open(Path(args.datasets), "r") as file:
        datasets = json.load(file)

    taxon_id = str(args.taxon_id)

    # 1) Try direct lookup in the local JSON first (no network)
    clade_match = None
    if isinstance(datasets, dict) and taxon_id in datasets:
        clade_match = datasets[taxon_id]

    # 2) If not found, fall back to NCBI lineage lookup
    if not clade_match:
        ncbi_url = f"{args.ncbi_url}/{taxon_id}/dataset_report"
        clade_match = get_dataset_match(ncbi_url, datasets)

    if not clade_match:
        # At this point, either the taxon really isn't covered,
        # or NCBI/network failed and we couldn't walk the lineage.
        raise ValueError("No match found")

    if args.output == "stdout":
        print(clade_match)
    else:
        with open(args.output, "w+") as output:
            # NOTE: this branch looks a bit odd, but I'm keeping your logic.
            if clade_match[0] == args.species:
                output.write(clade_match[1])
            else:
                output.write(clade_match[0])

    return None



if __name__ == "__main__":
    main()
