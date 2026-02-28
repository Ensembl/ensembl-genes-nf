#!/usr/bin/env python3
import duckdb, hashlib
from pathlib import Path

def connect(db_path: str):
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(db_path)
    # Minimal schema to start
    con.execute("CREATE TABLE IF NOT EXISTS metrics_core (run_id TEXT, sample_id TEXT, module TEXT, metric TEXT, scope TEXT, length INTEGER, value_num DOUBLE, unit TEXT)")
    con.execute("CREATE TABLE IF NOT EXISTS artifact (artifact_id TEXT, run_id TEXT, sample_id TEXT, module TEXT, name TEXT, path TEXT, checksum TEXT, size_bytes BIGINT, content_type TEXT)")
    con.execute("CREATE TABLE IF NOT EXISTS qc_rule_set (rule_set_id TEXT, name TEXT, version INTEGER, yaml_text TEXT, created_at TIMESTAMP)")
    con.execute("CREATE TABLE IF NOT EXISTS qc_eval (eval_id TEXT, run_id TEXT, sample_id TEXT, rule_set_id TEXT, scope TEXT, length INTEGER, rule_name TEXT, metric TEXT, op TEXT, threshold_num_min DOUBLE, threshold_num_max DOUBLE, observed_num DOUBLE, observed_text TEXT, severity TEXT, passed BOOLEAN, decided_at TIMESTAMP)")
    con.execute("CREATE TABLE IF NOT EXISTS gate_selection (run_id TEXT, sample_id TEXT, rule_set_id TEXT, selected_for_translon BOOLEAN, selected_for_trackhub BOOLEAN, selected_reason TEXT, pass_lengths_json TEXT)")
    return con

def upsert_rule_set(con, name: str, version: int, yaml_text: str) -> str:
    rule_set_id = f"{name}_v{version}"
    con.execute("INSERT INTO qc_rule_set(rule_set_id, name, version, yaml_text, created_at) VALUES (?, ?, ?, ?, current_timestamp)", [rule_set_id, name, version, yaml_text])
    return rule_set_id

def insert_metrics(con, rows):
    if not rows:
        return
    con.executemany("INSERT INTO metrics_core(run_id, sample_id, module, metric, scope, length, value_num, unit) VALUES (?, ?, ?, ?, ?, ?, ?, ?)", rows)

def insert_artifact(con, row):
    con.execute("INSERT INTO artifact(artifact_id, run_id, sample_id, module, name, path, checksum, size_bytes, content_type) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", row)

def md5sum(path):
    h = hashlib.md5()
    with open(path, rb) as f:
        for chunk in iter(lambda: f.read(8192), b):
            h.update(chunk)
    return h.hexdigest()
