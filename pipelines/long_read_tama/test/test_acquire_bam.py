import importlib.util
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
