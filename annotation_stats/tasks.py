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


def get_new_annotations(query_file_path):
    """
    Retrieve new annotation from the production metadata database.
    """
    with open(query_file_path, "r") as query_file:
        query = query_file.read()

    meta1 = {
        "host": "mysql-ens-meta-prod-1",
        "port": 4483,
        "user": "ensro",
        "database": "ensembl_metadata_qrp",
    }

    connection_args = meta1

    connection = pymysql.connect(
        host=connection_args["host"],
        port=connection_args["port"],
        user=connection_args["user"],
        database=connection_args["database"],
        # cursorclass=pymysql.cursors.DictCursor,
    )

    with connection:
        with connection.cursor() as cursor:
            cursor.execute(query)
            query_result = cursor.fetchall()

            columns = [column[0] for column in cursor.description]

    print(tabulate(query_result, headers=columns, tablefmt="psql"))

    exit()

    return query_result


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
