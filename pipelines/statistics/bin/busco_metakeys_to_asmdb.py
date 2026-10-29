#!/usr/bin/env python
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
"""Load BUSCO genome results into the assembly metadata database (assembly_metrics table)."""

import argparse
import json
import re
from pathlib import Path
from typing import Any, Dict, Optional, Union

import pymysql

GENOME_MODES = ("genome", "euk_genome_met", "euk_genome_min")
BUSCO_STATUS_EVENT = "genome_busco.status"


def parse_connection_params(json_string: str) -> Dict[str, Any]:
    """Parse the assembly metadata DB connection params from a JSON string."""
    try:
        params = json.loads(json_string)
    except json.JSONDecodeError as err:
        raise ValueError(f"--asm-metadata is not valid JSON: {err}") from err
    if not isinstance(params, dict):
        raise ValueError("--asm-metadata must be a JSON object")
    if "port" in params:
        params["port"] = int(params["port"])
    return params


def get_assembly_id(gca: str, metadata_params: Dict[str, Any]) -> int:
    """Get the assembly_id for a versioned GCA accession, required to insert into assembly_metrics."""
    query = """
        SELECT assembly_id
        FROM assembly
        WHERE CONCAT(gca_chain, '.', gca_version) = %s
    """
    conn = pymysql.connect(**metadata_params)
    try:
        with conn.cursor() as cursor:
            cursor.execute(query, (gca,))
            row = cursor.fetchone()
    finally:
        conn.close()

    if not row:
        raise ValueError(f"No assembly_id found for GCA: {gca}")
    print(f"Assembly id for {gca} is {row[0]}")
    return row[0]


def parse_busco_file(file_path: str) -> Dict[str, Union[str, int]]:
    """
    Parse a BUSCO genome-mode short summary file into assembly.* metrics.

    Args:
        file_path: Path to the BUSCO short summary file.

    Returns:
        Dictionary of metric name to value.

    Raises:
        ValueError: if the file is not a genome-mode summary or has no score line.
    """
    content = Path(file_path).read_text(encoding="utf-8")

    def _search(pattern: str) -> Optional[str]:
        match = re.search(pattern, content)
        return match.group(1) if match else None

    mode = _search(r"BUSCO was run in mode: ([\w_]+)")
    if mode not in GENOME_MODES:
        raise ValueError(f"Expected a BUSCO genome mode summary, got mode '{mode}' in {file_path}")

    # The erroneous (E) value is only reported by some BUSCO versions/modes
    score_pattern = (
        r"C:(\d+\.\d+)%\[S:(\d+\.\d+)%.*,D:(\d+\.\d+)%\],F:(\d+\.\d+)%.*,M:(\d+\.\d+)%,n:(\d+)"
        r"(?:,E:(\d+\.\d+)%)?"
    )
    score_match = re.search(score_pattern, content)
    if not score_match:
        raise ValueError(f"No BUSCO score line found in {file_path}")

    data: Dict[str, Union[str, int]] = {
        "assembly.busco_version": str(_search(r"BUSCO version is: (\d+\.\d+\.\d+)")),
        "assembly.busco_dataset": str(_search(r"The lineage dataset is: ([\w_]+)")),
        "assembly.busco": score_match.group(0),
        "assembly.busco_mode": "genome",
        "assembly.busco_completeness": str(_search(r"(\d+)\s+Complete BUSCOs \(C\)")),
        "assembly.busco_single_copy": str(_search(r"(\d+)\s+Complete and single-copy BUSCOs \(S\)")),
        "assembly.busco_duplicated": str(_search(r"(\d+)\s+Complete and duplicated BUSCOs \(D\)")),
        "assembly.busco_fragmented": str(_search(r"(\d+)\s+Fragmented BUSCOs \(F\)")),
        "assembly.busco_missing": str(_search(r"(\d+)\s+Missing BUSCOs \(M\)")),
        "assembly.busco_total": int(score_match.group(6)),
    }
    if score_match.group(7) is not None:
        data["assembly.busco_erroneus"] = score_match.group(7)
    return data


def insert_metrics(
    metrics: Dict[str, Union[str, int]], assembly_id: int, metadata_params: Dict[str, Any]
) -> None:
    """Insert or update BUSCO metrics in the assembly_metrics table for a given assembly_id."""
    query = """
        INSERT INTO assembly_metrics (assembly_id, metrics_name, metrics_value)
        VALUES (%s, %s, %s)
        ON DUPLICATE KEY UPDATE metrics_value = VALUES(metrics_value)
    """
    rows = [(assembly_id, key, value) for key, value in metrics.items()]
    conn = pymysql.connect(**metadata_params)
    try:
        with conn.cursor() as cursor:
            cursor.executemany(query, rows)
        conn.commit()
    finally:
        conn.close()
    print(f"Inserted {len(rows)} BUSCO metrics for assembly_id {assembly_id}")


def mark_busco_done(assembly_id: int, metadata_params: Dict[str, Any]) -> None:
    """Set genome_busco.status to done in the assembly_events table for a given assembly_id.

    Upserts, so an assembly with no genome_busco.status event yet (e.g. a manual run) gets one.
    """
    query = """
        INSERT INTO assembly_events (assembly_id, event, status)
        VALUES (%s, %s, %s)
        ON DUPLICATE KEY UPDATE status = VALUES(status)
    """
    conn = pymysql.connect(**metadata_params)
    try:
        with conn.cursor() as cursor:
            cursor.execute(query, (assembly_id, BUSCO_STATUS_EVENT, "done"))
        conn.commit()
    finally:
        conn.close()
    print(f"Marked {BUSCO_STATUS_EVENT} as done for assembly_id {assembly_id}")


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Parse a BUSCO genome summary and load it into the assembly metadata DB."
    )
    parser.add_argument("--version", action="version", version="busco_metakeys_to_asmdb.py version 1.0.0")
    parser.add_argument("--file", required=True, help="Path to the BUSCO short summary file")
    parser.add_argument("--gca", required=True, help="Versioned GCA accession")
    parser.add_argument(
        "--asm-metadata",
        required=True,
        help="JSON string with pymysql connection params (host, port, user, password, database)",
    )
    parser.add_argument("--output-json", help="Optional path to also write the parsed metrics as JSON")
    return parser


def main() -> None:
    """Main entry-point."""
    args = _build_arg_parser().parse_args()

    metadata_params = parse_connection_params(args.asm_metadata)
    metrics = parse_busco_file(args.file)

    if args.output_json:
        Path(args.output_json).write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    assembly_id = get_assembly_id(args.gca, metadata_params)
    insert_metrics(metrics, assembly_id, metadata_params)
    mark_busco_done(assembly_id, metadata_params)


if __name__ == "__main__":
    main()
