import importlib.util
import sys
from pathlib import Path


SPEC = importlib.util.spec_from_file_location(
    "download_fastq",
    Path(__file__).parents[1] / "bin" / "download_fastq.py",
)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def test_install_cache_file_handles_separate_filesystems(tmp_path):
    fetched = tmp_path / "work" / "reads.fastq.gz"
    cache = tmp_path / "cache" / "reads.fastq.gz"
    fetched.parent.mkdir()
    cache.parent.mkdir()
    fetched.write_bytes(b"downloaded-fastq")

    mod.install_cache_file(fetched, cache)

    assert cache.read_bytes() == b"downloaded-fastq"
    assert not fetched.exists()
    assert not Path(str(cache) + ".partial").exists()


def test_link_cached_file_does_not_copy_contents(tmp_path):
    cache = tmp_path / "cache" / "reads.fastq.gz"
    output = tmp_path / "work" / "downloaded" / cache.name
    cache.parent.mkdir()
    cache.write_bytes(b"cached-fastq")

    mod.link_cached_file(cache, output)

    assert output.is_symlink()
    assert output.resolve() == cache.resolve()
    assert output.read_bytes() == cache.read_bytes()


def test_cache_hit_links_without_rehashing(tmp_path, monkeypatch):
    expected = "0123456789abcdef0123456789abcdef"
    cache = tmp_path / "cache" / "ONT_FASTQ" / expected / "reads.fastq.gz"
    cache.parent.mkdir(parents=True)
    cache.write_bytes(b"cached-fastq")
    output = tmp_path / "work" / "downloaded"
    report = output / "checksum.tsv"

    monkeypatch.setattr(mod, "validate_gzip", lambda _: (_ for _ in ()).throw(AssertionError("cache revalidated")))
    monkeypatch.setattr(mod, "checksum", lambda _: (_ for _ in ()).throw(AssertionError("cache rehashed")))
    monkeypatch.setattr(sys, "argv", [
        "download_fastq.py", "RUN", expected, cache.name, "", str(tmp_path / "cache"),
        "ONT_FASTQ", str(output), "1", str(report),
    ])

    mod.main()

    assert (output / cache.name).is_symlink()
    assert (output / cache.name).resolve() == cache.resolve()
