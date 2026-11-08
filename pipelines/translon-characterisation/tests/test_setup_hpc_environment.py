import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / 'bin' / 'setup_hpc_environment.py'
REPO = Path(__file__).parents[3]


class HpcSetupTest(unittest.TestCase):
    def test_setup_is_idempotent_and_writes_explicit_config(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / 'translon'
            result = subprocess.run([
                sys.executable, str(SCRIPT), '--root', str(root), '--repo', str(REPO),
                '--partition', 'standard', '--account', 'project-a', '--skip-runtime-check',
            ], check=True, capture_output=True, text=True)
            self.assertIn('Prepared', result.stdout)
            for name in ('apptainer-cache', 'tmp', 'work', 'logs', 'results', 'inputs', 'references', 'models'):
                self.assertTrue((root / name).is_dir(), name)
            config = (root / 'hpc.config').read_text()
            self.assertIn("queue = 'standard'", config)
            self.assertIn("clusterOptions = '--account=project-a'", config)
            self.assertIn("cacheDir = '/", config)
            self.assertFalse((root / 'hpc.env').exists())
            self.assertFalse((root / 'run_characterisation.sh').exists())
            report = json.loads((root / 'setup-report.json').read_text())
            self.assertEqual(report['profile'], 'slurm_apptainer')
            self.assertEqual(report['config'], str((root / 'hpc.config').resolve()))

            second = subprocess.run([
                sys.executable, str(SCRIPT), '--root', str(root), '--repo', str(REPO),
                '--skip-runtime-check',
            ], check=True, capture_output=True, text=True)
            self.assertIn('Prepared', second.stdout)

    def test_required_input_manifest_is_checked(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / 'translon'
            manifest = Path(tmp) / 'inputs.tsv'
            manifest.write_text('name\tpath\trequired\nintervals\t' + str(Path(tmp) / 'missing.bed') + '\tyes\n')
            result = subprocess.run([
                sys.executable, str(SCRIPT), '--root', str(root), '--repo', str(REPO),
                '--inputs-manifest', str(manifest), '--skip-runtime-check',
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertIn('required input is missing', result.stderr)


if __name__ == '__main__':
    unittest.main()
