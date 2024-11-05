#!/usr/bin/env python
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""A simple helper script to connect to a core database and retrieve a single meta_value 
or multiple values dumped to json."""

__all__ = ["get_single_meta_value", "get_multi_meta_values"]

from pathlib import Path
import logging

from sqlalchemy.engine import URL

from ensembl.io.genomio.utils.json_utils import print_json
from ensembl.io.genomio.database import DBConnectionLite
from ensembl.utils.argparse import ArgumentParser
from ensembl.utils.logging import init_logging_with_args


def get_single_meta_value(server_url: URL, db_name: str, single_meta_key: str) -> str:
    """Returns a single meta_value from a core database. The first meta_value encountered on
    this meta_key will be returned.

    Args:
        server_url: Server URL where the core databases are stored.
        db_name: Name of the target DB to query.
        single_meta_key: The meta table 'meta_key' to query.

    Returns:
        meta_value: Value obtained from the meta table meta_key pair.

    """
    db_url = server_url.set(database=db_name)
    core_db = DBConnectionLite(db_url)
    meta_value = core_db.get_meta_value(f"{single_meta_key}")

    return meta_value


def get_multi_meta_values(server_url: URL, db_name: str, query_meta_keys: Path) -> dict[str, str]:
    """Returns a set of values based on set of 2 or more input DB meta_keys.

    Args:
        server_url: Server URL where the core databases are stored.
        db_name: Name of the target DB to query.
        query_meta_keys: The meta table 'meta_key' list to query.

    """
    db_url = server_url.set(database=db_name)
    core_db = DBConnectionLite(db_url)
    query_meta_values = {}
    input_metakey_count = 0
    unpopulated_meta_keys = []

    with Path(query_meta_keys).open("r") as fh:
        for file_inline in fh:
            input_metakey_count += 1
            meta_key = file_inline.strip()
            meta_value = core_db.get_meta_value(f"{meta_key}")

            if meta_value is not None:
                query_meta_values[f"{meta_key}"] = meta_value
            else:
                unpopulated_meta_keys.append(f"{meta_key}")
                logging.info(f"Meta query returned no entry on meta_key: '{meta_key}'")

    # Now assess what meta info was recovered and dump to JSON
    total_queries_located = len(query_meta_values.items())
    meta_populated = False

    if total_queries_located == input_metakey_count:
        meta_populated = True
    elif (total_queries_located >= 1) and (total_queries_located < input_metakey_count):
        meta_populated = True
        logging.info(
            f"Not all query meta_keys were located [Queries (input: {input_metakey_count}) vs (Located: {total_queries_located})"
        )
        logging.info(f"Meta_key query / queries missing: {unpopulated_meta_keys}")
    else:
        logging.warning("None of the input query meta keys were located. No data to dump to JSON")

    if meta_populated is True:
        query_meta_values["database_name"] = f"{db_name}"
        print_json("coredb_meta.json", query_meta_values)

    return


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_server_arguments()
    parser.add_argument("--database_name", default=None, help="Target database name.")
    parser.add_argument("--meta_key", default=None, help="Single meta_key to query target database.")
    parser.add_argument(
        "--meta_keys_list", default=None, help="Input file with >=2 meta_keys to query target database."
    )
    parser.add_log_arguments(add_log_file=False)
    args = parser.parse_args()
    init_logging_with_args(args)

    get_multi_meta_values(args.url, args.database_name, args.meta_keys_list)


if __name__ == "__main__":
    main()
