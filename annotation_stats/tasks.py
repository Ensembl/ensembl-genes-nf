#!/usr/bin/env python3
# -*- coding: utf-8 -*-


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
    fire.Fire()
