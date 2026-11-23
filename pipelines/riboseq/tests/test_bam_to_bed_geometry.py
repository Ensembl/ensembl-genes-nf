import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "bin" / "bam_to_bed.py"
SPEC = importlib.util.spec_from_file_location("bam_to_bed", SCRIPT)
bam_to_bed = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bam_to_bed)


def test_forward_reverse_and_five_prime_clips():
    length, clip, positions = bam_to_bed.cigar_geometry([(4, 3), (0, 28)], 100, False)
    assert (length, clip, positions[0], positions[-1]) == (28, 3, 100, 127)
    length, clip, positions = bam_to_bed.cigar_geometry([(0, 28), (4, 4)], 100, True)
    assert (length, clip, positions[0], positions[-1]) == (28, 4, 127, 100)
    assert bam_to_bed.is_profile_eligible(20, 0)
    assert bam_to_bed.is_profile_eligible(24, 0)
    assert not bam_to_bed.is_profile_eligible(20, 1)


def test_insertions_deletions_and_splices_project_query_offsets():
    _, _, insertion = bam_to_bed.cigar_geometry([(0, 10), (1, 2), (0, 10)], 100, False)
    assert insertion[10] is None
    assert insertion[12] == 110
    _, _, deletion = bam_to_bed.cigar_geometry([(0, 10), (2, 2), (0, 10)], 100, False)
    assert deletion[10] == 112
    _, _, splice = bam_to_bed.cigar_geometry([(0, 10), (3, 100), (0, 10)], 100, False)
    assert splice[10] == 210
