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
Task functions for the annotation stats Nextflow pipeline.
"""


# standard library
import pathlib

# third party
import fire

# project


def check_stats_files(annotation_directory: str, production_name: str):
    annotation_directory = pathlib.Path(annotation_directory)

    readme_file = annotation_directory / "statistics_README.txt"
    statistics_file = annotation_directory / f"{production_name}_annotation_statistics.txt"

    return readme_file.exists() and statistics_file.exists()


if __name__ == "__main__":
    try:
        fire.Fire()
    except KeyboardInterrupt:
        print("Interrupted with CTRL-C, exiting...")
