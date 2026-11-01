#!/usr/bin/env python3
"""Write the small, fixed-schema evidence relation table as Parquet.

This is deliberately limited to the relation table emitted by the long-read
evidence packager: required UTF-8 strings, one row group, uncompressed PLAIN
pages, and no optional statistics.  Keeping this writer dependency-free makes
the evidence process reproducible on hosts where pyarrow/polars/duckdb are not
installed while still producing a genuine Parquet file.
"""
from __future__ import annotations

import struct
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence


PARQUET_MAGIC = b"PAR1"
BYTE_ARRAY = 6
REQUIRED = 0
PLAIN = 0
RLE = 3
UNCOMPRESSED = 0
DATA_PAGE = 0

# TCompactProtocol type identifiers.
STOP = 0
BOOLEAN_TRUE = 1
I32 = 5
I64 = 6
BINARY = 8
LIST = 9
STRUCT = 12


class CompactWriter:
    """Minimal Apache Thrift compact-protocol writer used by Parquet metadata."""

    def __init__(self) -> None:
        self.data = bytearray()
        self._last_field_ids: list[int] = []
        self._last_field_id = 0

    def byte(self, value: int) -> None:
        self.data.append(value & 0xFF)

    def raw(self, value: bytes) -> None:
        self.data.extend(value)

    def varint(self, value: int) -> None:
        if value < 0:
            raise ValueError("compact-protocol varints must be non-negative")
        while value & ~0x7F:
            self.byte((value & 0x7F) | 0x80)
            value >>= 7
        self.byte(value)

    def zigzag32(self, value: int) -> int:
        return (value << 1) ^ (value >> 31)

    def zigzag64(self, value: int) -> int:
        return (value << 1) ^ (value >> 63)

    def i32(self, value: int) -> None:
        self.varint(self.zigzag32(value))

    def i64(self, value: int) -> None:
        self.varint(self.zigzag64(value))

    def binary(self, value: bytes) -> None:
        self.varint(len(value))
        self.raw(value)

    def field(self, field_id: int, field_type: int) -> None:
        delta = field_id - self._last_field_id
        if 1 <= delta <= 15:
            self.byte((delta << 4) | field_type)
        else:
            self.byte(field_type)
            self.i32(field_id)
        self._last_field_id = field_id

    def begin_struct(self) -> None:
        self._last_field_ids.append(self._last_field_id)
        self._last_field_id = 0

    def end_struct(self) -> None:
        self.byte(STOP)
        self._last_field_id = self._last_field_ids.pop()

    def list(self, element_type: int, size: int, write_element: Callable[[], None]) -> None:
        if size < 15:
            self.byte((size << 4) | element_type)
        else:
            self.byte(0xF0 | element_type)
            self.varint(size)
        for _ in range(size):
            write_element()


def _field(writer: CompactWriter, field_id: int, field_type: int, write_value: Callable[[], None]) -> None:
    writer.field(field_id, field_type)
    write_value()


def _string(writer: CompactWriter, value: str) -> None:
    writer.binary(value.encode("utf-8"))


def _schema_element(writer: CompactWriter, name: str, *, root: bool, child_count: int = 0) -> None:
    writer.begin_struct()
    if not root:
        _field(writer, 1, I32, lambda: writer.i32(BYTE_ARRAY))
        _field(writer, 3, I32, lambda: writer.i32(REQUIRED))
    _field(writer, 4, BINARY, lambda: _string(writer, name))
    if root:
        _field(writer, 5, I32, lambda: writer.i32(child_count))
    else:
        # ConvertedType.UTF8 is enum value 0.
        _field(writer, 6, I32, lambda: writer.i32(0))
    writer.end_struct()


def _data_page_header(num_values: int, uncompressed_size: int, minimum: str | None, maximum: str | None) -> bytes:
    writer = CompactWriter()
    writer.begin_struct()
    _field(writer, 1, I32, lambda: writer.i32(DATA_PAGE))
    _field(writer, 2, I32, lambda: writer.i32(uncompressed_size))
    _field(writer, 3, I32, lambda: writer.i32(uncompressed_size))
    _field(writer, 5, STRUCT, lambda: _write_data_page_header(writer, num_values, minimum, maximum))
    writer.end_struct()
    return bytes(writer.data)


def _write_data_page_header(writer: CompactWriter, num_values: int, minimum: str | None, maximum: str | None) -> None:
    writer.begin_struct()
    _field(writer, 1, I32, lambda: writer.i32(num_values))
    _field(writer, 2, I32, lambda: writer.i32(PLAIN))
    _field(writer, 3, I32, lambda: writer.i32(RLE))
    _field(writer, 4, I32, lambda: writer.i32(RLE))
    if minimum is not None and maximum is not None:
        _field(writer, 5, STRUCT, lambda: _statistics(writer, minimum, maximum))
    writer.end_struct()


def _statistics(writer: CompactWriter, minimum: str, maximum: str) -> None:
    writer.begin_struct()
    _field(writer, 1, BINARY, lambda: _string(writer, maximum))
    _field(writer, 2, BINARY, lambda: _string(writer, minimum))
    _field(writer, 3, I64, lambda: writer.i64(0))
    _field(writer, 5, BINARY, lambda: _string(writer, maximum))
    _field(writer, 6, BINARY, lambda: _string(writer, minimum))
    _field(writer, 7, BOOLEAN_TRUE, lambda: None)
    _field(writer, 8, BOOLEAN_TRUE, lambda: None)
    writer.end_struct()


def _column_metadata(
    writer: CompactWriter,
    name: str,
    row_count: int,
    page_size: int,
    page_offset: int,
    minimum: str | None,
    maximum: str | None,
) -> None:
    writer.begin_struct()
    _field(writer, 1, I32, lambda: writer.i32(BYTE_ARRAY))

    encodings = [PLAIN, RLE]

    def write_encoding() -> None:
        writer.i32(encodings[write_encoding.index])
        write_encoding.index += 1

    write_encoding.index = 0  # type: ignore[attr-defined]
    _field(writer, 2, LIST, lambda: writer.list(I32, len(encodings), write_encoding))
    _field(writer, 3, LIST, lambda: writer.list(BINARY, 1, lambda: _string(writer, name)))
    _field(writer, 4, I32, lambda: writer.i32(UNCOMPRESSED))
    _field(writer, 5, I64, lambda: writer.i64(row_count))
    _field(writer, 6, I64, lambda: writer.i64(page_size))
    _field(writer, 7, I64, lambda: writer.i64(page_size))
    _field(writer, 9, I64, lambda: writer.i64(page_offset))
    if minimum is not None and maximum is not None:
        _field(writer, 12, STRUCT, lambda: _statistics(writer, minimum, maximum))
    writer.end_struct()


def _column_chunk(
    writer: CompactWriter,
    name: str,
    row_count: int,
    page_size: int,
    page_offset: int,
    minimum: str | None,
    maximum: str | None,
) -> None:
    writer.begin_struct()
    _field(writer, 2, I64, lambda: writer.i64(page_offset))
    _field(writer, 3, STRUCT, lambda: _column_metadata(writer, name, row_count, page_size, page_offset, minimum, maximum))
    writer.end_struct()


def _row_group(writer: CompactWriter, columns: Sequence[tuple[str, int, int, int, str | None, str | None]], row_count: int) -> None:
    writer.begin_struct()

    def write_column() -> None:
        name, page_size, page_offset, _, minimum, maximum = columns[write_column.index]
        write_column.index += 1
        _column_chunk(writer, name, row_count, page_size, page_offset, minimum, maximum)

    write_column.index = 0  # type: ignore[attr-defined]
    _field(writer, 1, LIST, lambda: writer.list(STRUCT, len(columns), write_column))
    _field(writer, 2, I64, lambda: writer.i64(sum(page_size for _, page_size, _, _, _, _ in columns)))
    _field(writer, 3, I64, lambda: writer.i64(row_count))
    writer.end_struct()


def _file_metadata(
    writer: CompactWriter,
    fields: Sequence[str],
    row_count: int,
    columns: Sequence[tuple[str, int, int, int, str | None, str | None]],
) -> None:
    writer.begin_struct()
    _field(writer, 1, I32, lambda: writer.i32(1))

    schema_elements = [("schema", True)] + [(name, False) for name in fields]

    def write_schema() -> None:
        name, root = schema_elements[write_schema.index]
        write_schema.index += 1
        _schema_element(writer, name, root=root, child_count=len(fields))

    write_schema.index = 0  # type: ignore[attr-defined]

    _field(writer, 2, LIST, lambda: writer.list(STRUCT, len(fields) + 1, write_schema))
    _field(writer, 3, I64, lambda: writer.i64(row_count))
    _field(writer, 4, LIST, lambda: writer.list(STRUCT, 1, lambda: _row_group(writer, columns, row_count)))
    _field(writer, 6, BINARY, lambda: _string(writer, "ensembl-genes-nf long-read evidence"))
    writer.end_struct()


def _encode_column(values: Iterable[str]) -> bytes:
    output = bytearray()
    for value in values:
        encoded = value.encode("utf-8")
        output.extend(struct.pack("<I", len(encoded)))
        output.extend(encoded)
    return bytes(output)


def write_string_table(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
    """Write ``rows`` with the supplied fields to a valid Parquet file."""
    materialized = [{field: str(row[field]) for field in fields} for row in rows]
    path.parent.mkdir(parents=True, exist_ok=True)

    data = bytearray(PARQUET_MAGIC)
    columns: list[tuple[str, int, int, int, str | None, str | None]] = []
    for field in fields:
        values = [row[field] for row in materialized]
        encoded = _encode_column(values)
        minimum = min(values) if values else None
        maximum = max(values) if values else None
        page_header = _data_page_header(len(values), len(encoded), minimum, maximum)
        page_offset = len(data)
        data.extend(page_header)
        data.extend(encoded)
        columns.append((field, len(page_header) + len(encoded), page_offset, len(encoded), minimum, maximum))

    footer_writer = CompactWriter()
    _file_metadata(footer_writer, fields, len(materialized), columns)
    footer = bytes(footer_writer.data)
    data.extend(footer)
    data.extend(struct.pack("<I", len(footer)))
    data.extend(PARQUET_MAGIC)
    path.write_bytes(data)
