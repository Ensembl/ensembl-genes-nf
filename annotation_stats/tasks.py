#!/usr/bin/env python3
# -*- coding: utf-8 -*-


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


"""
Annotation stats Nextflow pipeline tasks.
"""


# standard library
import pathlib

from pprint import pp as pprint

# third party
import fire
import pymysql

from tabulate import tabulate

# project


meta1 = {
    "host": "mysql-ens-meta-prod-1",
    "port": 4483,
    "user": "ensro",
}


def run_sql_query(
    query_file: str,
    connection_config: dict = meta1,
    database: str = "ensembl_metadata_qrp",
    debug: bool = False,
):
    """
    Run the SQL query in query_file on the server and database defined in connection_config.
    """
    with open(query_file, "r") as file:
        query = file.read()

    connection = pymysql.connect(
        host=connection_config["host"],
        port=connection_config["port"],
        user=connection_config["user"],
        database=database,
        # cursorclass=pymysql.cursors.DictCursor,
    )

    with connection:
        with connection.cursor() as cursor:
            cursor.execute(query)
            query_result = cursor.fetchall()

            columns = [column[0] for column in cursor.description]

    if debug:
        print(tabulate(query_result, headers=columns, tablefmt="psql"))
        exit()

    return (columns, query_result)


def get_recent_annotations(query_file: str, annotations_csv: str):
    """
    Retrieve recent annotations from the production metadata database.
    """
    database = "ensembl_metadata_qrp"
    columns, query_result = run_sql_query(
        query_file=query_file, connection_config=meta1, database=database
    )

    annotation_databases = [annotation[0] for annotation in query_result]

    # save annotation databases list to a file
    with open(annotations_csv, "w") as file:
        for annotation_database in annotation_databases:
            file.write(f"{annotation_database}\n")

    # output annotation databases list to stdout
    for annotation_database in annotation_databases:
        print(annotation_database)


def check_stats_files(annotation_directory: str, production_name: str):
    annotation_directory = pathlib.Path(annotation_directory)

    readme_file = annotation_directory / "statistics_README.txt"
    statistics_file = (
        annotation_directory / f"{production_name}_annotation_statistics.txt"
    )

    return readme_file.exists() and statistics_file.exists()


if __name__ == "__main__":
    try:
        fire.Fire()
    except KeyboardInterrupt:
        print("Interrupted with CTRL-C, exiting...")
