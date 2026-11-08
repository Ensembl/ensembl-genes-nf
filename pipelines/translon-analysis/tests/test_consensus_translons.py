import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / 'bin' / 'consensus_translons.py'


class ConsensusContractTest(unittest.TestCase):
    def test_two_callers_create_trusted_interval_and_peptide(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = root / 'inputs'
            inputs.mkdir()
            (root / 'genome.fa').write_text('>chr1\nATGAAACCCGGGTAA\n')
            header = 'sample_id\ttool\tchrom\tstart\tend\tstrand\tframe\ttranscript_id\textra_json\n'
            row = 'S1\t{tool}\tchr1\t0\t15\t+\t0\ttx1\t' + json.dumps({'_blocks': [[0, 5], [7, 15]]}) + '\n'
            (inputs / 'ribocode.tsv').write_text(header + row.replace('{tool}', 'ribocode'))
            (inputs / 'ribotricer.tsv').write_text(header + row.replace('{tool}', 'ribotricer'))

            subprocess.run([
                sys.executable, str(SCRIPT), '--input-dir', str(inputs),
                '--genome-fasta', str(root / 'genome.fa'),
                '--min-caller-agreement', '2',
                '--candidates', str(root / 'candidates.tsv'),
                '--bed12', str(root / 'candidates.bed12'),
                '--intervals', str(root / 'intervals.tsv'),
                '--verdicts', str(root / 'verdicts.tsv'),
            ], check=True)

            with (root / 'candidates.tsv').open() as handle:
                candidate = next(csv.DictReader(handle, delimiter='\t'))
            with (root / 'verdicts.tsv').open() as handle:
                verdict = next(csv.DictReader(handle, delimiter='\t'))

            self.assertEqual(candidate['status'], 'trusted')
            self.assertEqual(candidate['caller_count'], '2')
            self.assertEqual(candidate['transcript_ids'], 'tx1')
            self.assertEqual(json.loads(candidate['candidate_blocks_json'])[0]['blocks'], [[0, 5], [7, 15]])
            self.assertEqual(verdict['peptide_sequence'], 'MKPG*')


if __name__ == '__main__':
    unittest.main()

