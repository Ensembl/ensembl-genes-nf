import importlib.util
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
