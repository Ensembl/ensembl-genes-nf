import csv
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
            header = 'sample_id\ttool\tchrom\tstart\tend\tstrand\tframe\n'
            row = 'S1\t{tool}\tchr1\t0\t15\t+\t0\n'
            (inputs / 'ribocode.tsv').write_text(header + row.format(tool='ribocode'))
            (inputs / 'ribotricer.tsv').write_text(header + row.format(tool='ribotricer'))

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
            self.assertEqual(verdict['peptide_sequence'], 'MKPG*')


if __name__ == '__main__':
    unittest.main()

