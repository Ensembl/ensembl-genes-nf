#!/usr/bin/env python3
"""Join a per-instance axis payload to the original attributed instances."""
import argparse, json
def read(path):
    with open(path) as handle: yield from (json.loads(line) for line in handle if line.strip())
def main():
    p=argparse.ArgumentParser(); p.add_argument('--instances',required=True); p.add_argument('--axis-records',required=True); p.add_argument('--output',required=True); a=p.parse_args()
    payloads={x['instance_id']:x['axis'] for x in read(a.axis_records)}
    with open(a.output,'w') as out:
        for instance in read(a.instances):
            instance['axes']=instance.get('axes',{}); instance['axes'].update(payloads.get(instance['instance_id'],{}))
            out.write(json.dumps(instance,sort_keys=True)+'\n')
if __name__=='__main__': main()
