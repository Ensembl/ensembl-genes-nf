import importlib.util
import sys
from pathlib import Path


SPEC = importlib.util.spec_from_file_location(
    "acquire_bam",
    Path(__file__).parents[1] / "bin" / "acquire_bam.py",
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_link_cached_file_does_not_copy_contents(tmp_path):
    cache = tmp_path / "cache" / "reads.bam"
    output = tmp_path / "work" / "downloaded" / cache.name
    cache.parent.mkdir()
    cache.write_bytes(b"cached-bam")

    MODULE.link_cached_file(cache, output)

    assert output.is_symlink()
    assert output.resolve() == cache.resolve()
    assert output.read_bytes() == cache.read_bytes()


def test_cache_hit_links_without_rehashing(tmp_path, monkeypatch):
    expected = "0123456789abcdef0123456789abcdef"
    cache = tmp_path / "cache" / "PACBIO_CCS_BAM" / expected / "reads.bam"
    cache.parent.mkdir(parents=True)
    cache.write_bytes(b"cached-bam")
    output = tmp_path / "work" / "downloaded"
    report = output / "checksum.tsv"

    monkeypatch.setattr(MODULE, "md5", lambda _: (_ for _ in ()).throw(AssertionError("cache rehashed")))
    monkeypatch.setattr(sys, "argv", [
        "acquire_bam.py", expected, cache.name, "", "", "",
        str(tmp_path / "cache"), "PACBIO_CCS_BAM", "RUN", str(output), str(report),
    ])

    MODULE.main()

    assert (output / cache.name).is_symlink()
    assert (output / cache.name).resolve() == cache.resolve()
