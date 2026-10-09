from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1] / "bin"),
)

import rewrite_target_gff3


def test_rewrites_ids_parents_and_versions(
    tmp_path: Path,
) -> None:
    id_map = tmp_path / "id_map.tsv"
    id_map.write_text(
        "old_id\tnew_id\tnew_version\n"
        "OLDG\tNEWG\t2\n"
        "OLDT\tNEWT\t3\n"
        "OLDE\tNEWE\t4\n"
        "OLDP\tNEWP\t5\n",
        encoding="utf-8",
    )

    input_gff = tmp_path / "target.gff3"
    input_gff.write_text(
        "##gff-version 3\n"
        "1\t.\tgene\t1\t100\t.\t+\t.\t"
        "ID=gene:OLDG;gene_id=OLDG;version=1\n"
        "1\t.\tmRNA\t1\t100\t.\t+\t.\t"
        "ID=transcript:OLDT;Parent=gene:OLDG;"
        "transcript_id=OLDT;version=1\n"
        "1\t.\texon\t1\t20\t.\t+\t.\t"
        "Parent=transcript:OLDT;Name=OLDE;exon_id=OLDE;version=1\n"
        "1\t.\tCDS\t1\t20\t.\t+\t0\t"
        "ID=CDS:OLDP;Parent=transcript:OLDT;"
        "protein_id=OLDP;version=1\n"
        "1\t.\tbiological_region\t101\t200\t.\t+\t.\t"
        "logic_name=cpg\n",
        encoding="utf-8",
    )

    output_gff = tmp_path / "rewritten.gff3"

    updates = rewrite_target_gff3.load_id_map(id_map)

    (
        feature_lines,
        id_replacements,
        version_replacements,
        used_ids,
    ) = rewrite_target_gff3.rewrite_gff(
        input_gff,
        output_gff,
        updates,
    )

    assert output_gff.read_text(encoding="utf-8") == (
        "##gff-version 3\n"
        "1\t.\tgene\t1\t100\t.\t+\t.\t"
        "ID=gene:NEWG;gene_id=NEWG;version=2\n"
        "1\t.\tmRNA\t1\t100\t.\t+\t.\t"
        "ID=transcript:NEWT;Parent=gene:NEWG;"
        "transcript_id=NEWT;version=3\n"
        "1\t.\texon\t1\t20\t.\t+\t.\t"
        "Parent=transcript:NEWT;Name=NEWE;exon_id=NEWE;version=4\n"
        "1\t.\tCDS\t1\t20\t.\t+\t0\t"
        "ID=CDS:NEWP;Parent=transcript:NEWT;"
        "protein_id=NEWP;version=5\n"
        "1\t.\tbiological_region\t101\t200\t.\t+\t.\t"
        "logic_name=cpg\n"
    )

    assert feature_lines == 5
    assert id_replacements == 11
    assert version_replacements == 4
    assert used_ids == {"OLDG", "OLDT", "OLDE", "OLDP"}

def test_load_id_map_accepts_mapping_decisions_tsv(
    tmp_path: Path,
) -> None:
    decisions_tsv = tmp_path / "decisions.tsv"
    decisions_tsv.write_text(
        "type\taction\tcurrent_stable_id\tnew_stable_id\tnew_version\n"
        "gene\tmapped\tTARGETG\tREFG\t2\n"
        "gene\tmissing\t\t\t0\n"
        "gene\tnew\tTARGET_NEW\tNEWG\t1\n",
        encoding="utf-8",
    )

    assert rewrite_target_gff3.load_id_map(decisions_tsv) == {
        "TARGETG": ("REFG", "2"),
        "TARGET_NEW": ("NEWG", "1"),
    }