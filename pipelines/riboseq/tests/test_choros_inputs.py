import csv
import importlib.util
from pathlib import Path
import subprocess

import pysam


SCRIPT = Path(__file__).parents[1] / "bin" / "prepare_choros_inputs.py"
SPEC = importlib.util.spec_from_file_location("prepare_choros_inputs", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_frame_offsets_place_asite_in_frame_zero():
    offsets = [MODULE.nearest_frame_offset(15, frame, 28) for frame in range(3)]
    assert offsets == [15, 14, 16]
    assert all((frame + offset) % 3 == 0 for frame, offset in enumerate(offsets))


def test_convert_ribometric_annotation(tmp_path):
    annotation = tmp_path / "annotation.tsv"
    annotation.write_text(
        "transcript_id\tcds_start\tcds_end\ttranscript_length\n"
        "tx1\t10\t40\t50\n"
    )
    output = tmp_path / "lengths.tsv"
    MODULE.write_choros_lengths(annotation, output)
    with output.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows == [{
        "transcript": "tx1",
        "utr5_length": "10",
        "cds_length": "30",
        "utr3_length": "10",
    }]


def test_write_choros_offsets(tmp_path):
    output = tmp_path / "offsets.tsv"
    MODULE.write_choros_offsets({28: 15}, output)
    assert output.read_text() == "length\tframe_0\tframe_1\tframe_2\n28\t15\t14\t16\n"


def test_prepare_bam_preserves_collapsed_weight_and_drops_ambiguous(tmp_path):
    source = tmp_path / "grouped.bam"
    header = {"HD": {"VN": "1.6", "SO": "queryname"}, "SQ": [{"SN": "tx1", "LN": 100}]}
    with pysam.AlignmentFile(source, "wb", header=header) as bam:
        for query_name, start in [("unique_x7", 5), ("ambiguous_x3", 10), ("ambiguous_x3", 20)]:
            read = pysam.AlignedSegment()
            read.query_name = query_name
            read.query_sequence = "A" * 28
            read.flag = 0
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 255
            read.cigar = ((0, 28),)
            read.query_qualities = pysam.qualitystring_to_array("I" * 28)
            bam.write(read)

    output = tmp_path / "choros.bam"
    metrics = tmp_path / "metrics.tsv"
    subprocess.run(
        [
            "python3", str(Path(__file__).parents[1] / "bin" / "prepare_choros_bam.py"),
            "--bam", str(source), "--output", str(output), "--metrics", str(metrics),
        ],
        check=True,
    )
    with pysam.AlignmentFile(output, "rb") as bam:
        retained = list(bam.fetch(until_eof=True))
    assert len(retained) == 1
    assert retained[0].query_name == "unique_x7"
    assert retained[0].get_tag("ZW") == 7.0
    assert "query_groups_ambiguous\t1" in metrics.read_text()
