"""Compatibility shims required by the pinned plastid Cython sources."""

import numpy as _numpy
import collections as _collections
import collections.abc as _collections_abc

if not hasattr(_numpy, "int"):
    _numpy.int = int
if not hasattr(_numpy, "long"):
    _numpy.long = int
if not hasattr(_numpy, "float"):
    _numpy.float = float
if not hasattr(_numpy, "bool"):
    _numpy.bool = bool

if not hasattr(_collections, "Iterable"):
    _collections.Iterable = _collections_abc.Iterable
