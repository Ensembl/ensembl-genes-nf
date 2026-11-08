import json
import hashlib
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).parents[1] / 'bin' / 'translon_characterise.py'
PHYLOCSF = Path(__file__).parents[1] / 'bin' / 'phylocsf_axis.py'
MSA_SETUP = Path(__file__).parents[1] / 'bin' / 'setup_maf_reference.py'
HISTORY = Path(__file__).parents[1] / 'bin' / 'evolutionary_history.py'


class CharacterisationContractTest(unittest.TestCase):
    def test_standard_bed12_is_accepted_without_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bed = root / 'orfs.bed'
            bed.write_text('chr1\t2\t8\torf-bed12\t0\t+\t2\t8\t0\t2\t3,3\t0,3\n')
            (root / 'gencode.gff3').write_text(
                '##gff-version 3\nchr1\tGENCODE\tmRNA\t1\t20\t.\t+\t.\tID=transcript:tx1\n'
                'chr1\tGENCODE\texon\t1\t10\t.\t+\t.\tID=exon:e1;Parent=transcript:tx1\n')
            (root / 'genome.fa').write_text('>chr1\nCCATGTAACCCCCCCCCCCC\n')
            subprocess.run([sys.executable, str(SCRIPT), 'substrate', '--intervals', str(bed),
                            '--gff', str(root / 'gencode.gff3'), '--genome', str(root / 'genome.fa'),
                            '--instances', str(root / 'instances.jsonl'), '--unhosted', str(root / 'unhosted.jsonl')], check=True)
            records = [json.loads(line) for line in (root / 'instances.jsonl').read_text().splitlines()]
            self.assertEqual(records[0]['interval_id'], 'orf-bed12')

    def test_spliced_context_and_unhosted_locus_are_distinguished(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'intervals.tsv').write_text(
                'interval_id\tchrom\tstart\tend\tstrand\tframe\ttranscript_ids\tcandidate_blocks_json\n'
                'orf1\tchr1\t2\t8\t+\t0\ttx1\t[{"transcript_id":"tx1","blocks":[[2,5],[5,8]]}]\n'
                'orphan\tchr9\t1\t4\t+\t0\t\t\n')
            (root / 'gencode.gff3').write_text(
                '##gff-version 3\nchr1\tGENCODE\tmRNA\t1\t20\t.\t+\t.\tID=transcript:tx1\n'
                'chr1\tGENCODE\texon\t1\t10\t.\t+\t.\tID=exon:e1;Parent=transcript:tx1\n'
                'chr1\tGENCODE\tCDS\t11\t20\t.\t+\t0\tID=CDS:c1;Parent=transcript:tx1\n')
            (root / 'genome.fa').write_text('>chr1\nCCATGTAACCCCCCCCCCCC\n')
            subprocess.run([sys.executable, str(SCRIPT), 'substrate', '--intervals', str(root / 'intervals.tsv'), '--gff', str(root / 'gencode.gff3'), '--genome', str(root / 'genome.fa'), '--instances', str(root / 'instances.jsonl'), '--unhosted', str(root / 'unhosted.jsonl')], check=True)
            instances = [json.loads(x) for x in (root / 'instances.jsonl').read_text().splitlines()]
            unhosted = [json.loads(x) for x in (root / 'unhosted.jsonl').read_text().splitlines()]
            self.assertEqual(instances[0]['class'], 'uORF')
            self.assertEqual(instances[0]['source_transcript_ids'], ['tx1'])
            self.assertEqual(instances[0]['source_candidate_blocks'][0]['transcript_id'], 'tx1')
            self.assertEqual(instances[0]['admissibility']['state'], 'positive')
            self.assertEqual(unhosted[0]['claim_hint'], 'translated-no-compatible-transcript')
            self.assertEqual(unhosted[0]['admissibility']['state'], 'unattributable')

    def test_missing_tool_result_is_uninformative_not_negative(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'input.jsonl').write_text(json.dumps({'instance_id': 'i1'}) + '\n')
            subprocess.run([sys.executable, str(SCRIPT), 'axis', '--axis', 'phylocsf', '--input', str(root / 'input.jsonl'), '--output', str(root / 'axis.jsonl')], check=True)
            result = json.loads((root / 'axis.jsonl').read_text())
            self.assertEqual(result['axes']['phylocsf']['state'], 'uninformative')

    def test_phylocsf_requires_identity_and_matched_null(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'alignment.fa').write_text('>hg38 instance_id=i1\nATGTAA\n>panTro6\nATGTAA\n')
            (root / 'raw.txt').write_text('2.5\n')
            (root / 'null.json').write_text(json.dumps({'positive_threshold': 1.0, 'provenance': 'test'}))
            subprocess.run([sys.executable, str(PHYLOCSF), '--alignment', str(root / 'alignment.fa'), '--raw', str(root / 'raw.txt'), '--matched-null', str(root / 'null.json'), '--output', str(root / 'axis.json')], check=True)
            result = json.loads((root / 'axis.json').read_text())
            self.assertEqual(result['instance_id'], 'i1')
            self.assertEqual(result['axis']['phylocsf']['state'], 'positive')

    def test_evolutionary_history_preserves_one_to_many_candidates_and_denominator(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'instances.jsonl').write_text(json.dumps({
                'instance_id': 'i1', 'interval_id': 'orf1', 'transcript_id': 'tx1',
                'class': 'uORF', 'chrom': 'chr1', 'start': 2, 'end': 8,
                'strand': '+', 'frame': '0',
            }) + '\n')
            (root / 'candidates.jsonl').write_text('\n'.join([
                json.dumps({'interval_id': 'orf1', 'species': 'mouse', 'candidate_gene': 'G1',
                            'candidate_transcript': 'T1', 'orthology_status': 'one_to_many',
                            'cds_status': 'aligned', 'leader_status': 'ambiguous',
                            'orf_status': 'unknown', 'translation_status': 'untested',
                            'chrom': 'chr1', 'start': 10, 'end': 20, 'strand': '+'}),
                json.dumps({'interval_id': 'orf1', 'species': 'mouse', 'candidate_gene': 'G2',
                            'candidate_transcript': 'T2', 'orthology_status': 'one_to_many',
                            'cds_status': 'partial', 'leader_status': 'missing',
                            'orf_status': 'not_observed', 'translation_status': 'untested',
                            'chrom': 'chr1', 'start': 20, 'end': 10, 'strand': '+'}),
            ]) + '\n')
            subprocess.run([sys.executable, str(HISTORY), '--instances', str(root / 'instances.jsonl'),
                            '--evidence', str(root / 'candidates.jsonl'), '--output', str(root / 'history.jsonl')], check=True)
            result = json.loads((root / 'history.jsonl').read_text())
            history = result['axes']['evolutionary_history']
            self.assertEqual(result['transcript_id'], 'tx1')
            self.assertEqual(history['state'], 'uninformative')
            self.assertEqual(history['candidate_state'], 'enumerated')
            self.assertEqual(history['summary']['candidate_count'], 2)
            self.assertEqual(history['species'][0]['candidate_count'], 2)
            self.assertEqual(history['summary']['species_states']['mouse']['gene'], 'ambiguous')
            self.assertEqual(history['summary']['species_states']['mouse']['transcript'], 'ambiguous')
            self.assertEqual(history['summary']['species_states']['mouse']['leader'], 'ambiguous')
            self.assertEqual(history['summary']['species_states']['mouse']['path'], 'ambiguous')
            self.assertFalse(history['denominator']['changed'])

    def test_phylocsf_batch_output_uses_explicit_identity_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'identities.jsonl').write_text(json.dumps({'alignment': 'one.msa.fasta', 'instance_id': 'i1'}) + '\n')
            (root / 'raw.txt').write_text('msas/one.msa.fasta\tdone\t2.5\n')
            (root / 'null.json').write_text(json.dumps({'positive_threshold': 1.0, 'provenance': 'test'}))
            subprocess.run([sys.executable, str(PHYLOCSF), '--identities', str(root / 'identities.jsonl'),
                            '--raw', str(root / 'raw.txt'), '--matched-null', str(root / 'null.json'),
                            '--output', str(root / 'axis.jsonl')], check=True)
            result = json.loads((root / 'axis.jsonl').read_text())
            self.assertEqual(result['instance_id'], 'i1')
            self.assertEqual(result['axis']['phylocsf']['state'], 'positive')

    def test_maf_setup_rejects_unchecked_reference_and_materialises_checked_pair(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            maf, index = root / 'chr1.maf.source', root / 'chr1.maf.bb.source'
            maf.write_text('##maf version=1\n')
            index.write_bytes(b'index')
            sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
            manifest = root / 'manifest.jsonl'
            manifest.write_text(json.dumps({'chrom': 'chr1', 'maf_url': maf.as_uri(), 'maf_sha256': sha(maf)}) + '\n')
            instances = root / 'instances.jsonl'
            instances.write_text(json.dumps({'chrom': 'chr1'}) + '\n')
            genome = root / 'genome.fa'
            genome.write_text('>chr1\n' + 'A' * 100 + '\n')
            fake_index = root / 'mafIndex'
            fake_index.write_text('#!/bin/sh\ntouch "$2"\n')
            fake_index.chmod(0o755)
            subprocess.run([sys.executable, str(MSA_SETUP), '--manifest', str(manifest),
                            '--instances', str(instances), '--genome', str(genome),
                            '--maf-index', str(fake_index), '--outdir', str(root / 'out'),
                            '--report', str(root / 'report.jsonl')], check=True)
            self.assertEqual((root / 'out' / 'chr1.maf').read_text(), maf.read_text())
            self.assertTrue((root / 'out' / 'chr1.maf.bb').is_file())


if __name__ == '__main__':
    unittest.main()
