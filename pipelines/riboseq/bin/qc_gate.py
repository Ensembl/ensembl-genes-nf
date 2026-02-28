#!/usr/bin/env python3
import argparse, json, yaml
from pathlib import Path
from qc_db import connect, upsert_rule_set

OPS = {
    '>=': lambda x, t: x is not None and x >= t,
    '<=': lambda x, t: x is not None and x <= t,
    'between': lambda x, t: x is not None and t[0] <= x <= t[1],
}

def get_metric(con, run_id, sample_id, metric, scope, length=None):
    if scope == 'sample':
        q = "SELECT value_num FROM metrics_core WHERE run_id=? AND sample_id=? AND metric=? AND scope='sample' LIMIT 1"
        r = con.execute(q, [run_id, sample_id, metric]).fetchone()
        return r[0] if r else None
    else:
        q = "SELECT value_num FROM metrics_core WHERE run_id=? AND sample_id=? AND metric=? AND scope='length' AND length=? LIMIT 1"
        r = con.execute(q, [run_id, sample_id, metric, int(length)]).fetchone()
        return r[0] if r else None

def read_offsets_table(path):
    rows = []
    with open(path) as f:
        header = f.readline()
        for line in f:
            if not line.strip():
                continue
            rl, off = line.split('	')[:2]
            rows.append((int(rl), int(off)))
    return rows

def write_filtered_offsets(in_path, out_path, keep_lengths):
    with open(in_path) as fi, open(out_path, 'w') as fo:
        header = fi.readline()
        fo.write(header)
        for line in fi:
            if not line.strip():
                continue
            L = int(line.split('	', 1)[0])
            if L in keep_lengths:
                fo.write(line)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--db', required=True)
    ap.add_argument('--run-id', required=True)
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--best-offset', required=True)
    ap.add_argument('--rules-yaml', required=True)
    ap.add_argument('--rule-set-name', default='default')
    ap.add_argument('--apply-for', default='translon,trackhub')
    ap.add_argument('--out-prefix', required=True)
    args = ap.parse_args()

    con = connect(args.db)
    rules_text = Path(args.rules_yaml).read_text()
    rules = yaml.safe_load(rules_text)
    rule_set_id = upsert_rule_set(con, rules.get('rule_set_name', args.rule_set_name), rules.get('version', 1), rules_text)

    offsets = read_offsets_table(args.best_offset)
    keep = set()
    eval_records = []

    for L, _O in offsets:
        passed_all = True
        for rule in rules.get('length_rules', []):
            metric = rule['metric']
            op = rule['op']
            thr = rule['threshold']
            val = get_metric(con, args.run_id, args.sample_id, metric, 'length', L)
            if val is None and not rule.get('include_if_missing', False):
                passed = False
            else:
                if isinstance(thr, (list, tuple)) and op == 'between':
                    passed = OPS['between'](val, thr)
                else:
                    passed = OPS[op](val, thr)
            eval_records.append((args.run_id, args.sample_id, rule_set_id, 'length', int(L), rule['name'], metric, op,
                                 float(thr) if isinstance(thr, (int, float)) else None,
                                 None,
                                 float(val) if val is not None else None,
                                 None, rule['severity'], bool(passed)))
            if not passed:
                passed_all = False
        if passed_all:
            keep.add(L)

    sample_pass = True
    for rule in rules.get('sample_rules', []):
        metric = rule['metric']
        op = rule['op']
        thr = rule['threshold']
        val = get_metric(con, args.run_id, args.sample_id, metric, 'sample')
        if val is None and not rule.get('include_if_missing', False):
            passed = False
        else:
            if isinstance(thr, (list, tuple)) and op == 'between':
                passed = OPS['between'](val, thr)
            else:
                passed = OPS[op](val, thr)
        eval_records.append((args.run_id, args.sample_id, rule_set_id, 'sample', None, rule['name'], metric, op,
                             float(thr) if isinstance(thr, (int, float)) else None, None,
                             float(val) if val is not None else None, None, rule['severity'], bool(passed)))
        if not passed:
            sample_pass = False

    con.executemany("""
        INSERT INTO qc_eval(eval_id, run_id, sample_id, rule_set_id, scope, length, rule_name, metric, op,
                            threshold_num_min, threshold_num_max, observed_num, observed_text, severity, passed, decided_at)
        VALUES (uuid(), ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, current_timestamp)
    """, eval_records)

    sel_t = 'translon' in args.apply_for and sample_pass and len(keep) > 0
    sel_h = 'trackhub' in args.apply_for and sample_pass and len(keep) > 0
    reason = 'sample rules pass and ≥1 passing length' if (sample_pass and keep) else 'rules not met'
    con.execute("""
        INSERT INTO gate_selection(run_id, sample_id, rule_set_id, selected_for_translon, selected_for_trackhub, selected_reason, pass_lengths_json)
        VALUES (?, ?, ?, ?, ?, ?, ?)
    """, [args.run_id, args.sample_id, rule_set_id, sel_t, sel_h, reason, json.dumps(sorted(list(keep)))])

    pass_lengths_path = f"{args.out_prefix}.pass_lengths.tsv"
    with open(pass_lengths_path, 'w') as f:
        f.write('length
')
        for L in sorted(keep):
            f.write(f'{L}
')

    filtered_offsets_path = f"{args.out_prefix}.offsets.pass.tsv"
    write_filtered_offsets(args.best_offset, filtered_offsets_path, keep)

    with open(f"{args.out_prefix}.qc.json", 'w') as f:
        json.dump({'sample_id': args.sample_id,
                   'rule_set_id': rule_set_id,
                   'selected_for_translon': sel_t,
                   'selected_for_trackhub': sel_h,
                   'pass_lengths': sorted(list(keep))}, f)

if __name__ == '__main__':
    main()
