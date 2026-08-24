# ORF-RATER Python 3 port

This directory is a maintained Python 3 port of the upstream ORF-RATER
sources at commit `5ae4cd69fc845e811e9c978ddbe48b706a071061`.

The port preserves the original command-line scripts and data formats. Changes
are limited to Python 3 language/runtime compatibility and removal of APIs
removed from current NumPy, pandas, pysam, and related dependencies. Numerical
parity must be established with a legacy Python 2 reference run before this
runtime replaces the legacy image in production.

The `multiisotonic/` package is vendored from upstream commit
`107a27870670212a6d981c05f8858fc33c8b0634` because that repository does not
contain install metadata.

plastid is built from pinned source commit
`d97f239d73b3a7c2eff46f71928b777431891f90` with its Cython extensions
enabled. The image applies compatibility fixes for removed NumPy aliases,
modern Cython types, and current C compiler headers; the complete reader
extension set is retained because plastid imports it at package startup.

`sitecustomize.py` supplies the legacy NumPy and `collections.Iterable` names
required by the pinned plastid runtime. This shim is isolated in the image and
does not change ORF-RATER's data formats.
