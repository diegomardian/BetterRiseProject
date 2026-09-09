"""Read native TIFF header facts without downloading an H&E image.

Only the beginning of a TIFF is needed to find its first image-file directory
and its pixel-size tags in the usual baseline-TIFF layout.  This module makes
no claim that a deposited image is the scanner's original full-resolution
object; it reports only what its native TIFF header declares.
"""

from __future__ import annotations

import struct
from collections.abc import Callable, Mapping


class TiffHeaderError(ValueError):
    """The bounded header prefix cannot establish TIFF dimensions or scale."""


_TYPE_FORMATS = {
    3: "H",  # SHORT
    4: "I",  # LONG
    5: "II",  # RATIONAL
    16: "Q",  # LONG8
}
_TYPE_SIZES = {3: 2, 4: 4, 5: 8, 16: 8}
_PIXEL_WIDTH = 256
_PIXEL_HEIGHT = 257
_X_RESOLUTION = 282
_Y_RESOLUTION = 283
_RESOLUTION_UNIT = 296
_REQUIRED_TAGS = {_PIXEL_WIDTH, _PIXEL_HEIGHT, _X_RESOLUTION, _Y_RESOLUTION, _RESOLUTION_UNIT}
_MAX_DIRECTORY_BYTES = 65_536


def _decode_value(
    raw: bytes, *, field_type: int, count: int, order: str
) -> int | tuple[int, int]:
    if count != 1 or field_type not in _TYPE_SIZES:
        raise TiffHeaderError(f"unsupported TIFF tag type/count: type={field_type}, count={count}")
    size = _TYPE_SIZES[field_type]
    if len(raw) != size:
        raise TiffHeaderError("TIFF range response was shorter than the requested tag value")
    return struct.unpack(f"{order}{_TYPE_FORMATS[field_type]}", raw)


def _header_layout(initial: bytes) -> tuple[str, int, int, int, int, int, str, str, str]:
    """Return byte order and first-directory layout from the initial TIFF bytes."""
    if len(initial) < 8:
        raise TiffHeaderError("header prefix is shorter than a TIFF header")
    endian = initial[:2]
    if endian == b"II":
        order = "<"
    elif endian == b"MM":
        order = ">"
    else:
        raise TiffHeaderError("not a TIFF byte-order marker")
    magic = struct.unpack_from(f"{order}H", initial, 2)[0]
    if magic == 42:
        ifd_offset = struct.unpack_from(f"{order}I", initial, 4)[0]
        count_size, count_field_size, entry_size, inline_size = 2, 4, 12, 4
        count_format, offset_format = "H", "I"
        kind = "classic_tiff"
    elif magic == 43:
        if len(initial) < 16:
            raise TiffHeaderError("header prefix is shorter than a BigTIFF header")
        offset_size, reserved = struct.unpack_from(f"{order}HH", initial, 4)
        if offset_size != 8 or reserved != 0:
            raise TiffHeaderError("unsupported BigTIFF offset layout")
        ifd_offset = struct.unpack_from(f"{order}Q", initial, 8)[0]
        count_size, count_field_size, entry_size, inline_size = 8, 8, 20, 8
        count_format, offset_format = "Q", "Q"
        kind = "big_tiff"
    else:
        raise TiffHeaderError(f"unsupported TIFF magic value {magic}")
    return (
        order,
        ifd_offset,
        count_size,
        count_field_size,
        entry_size,
        inline_size,
        count_format,
        offset_format,
        kind,
    )


def parse_tiff_header_ranges(
    read_range: Callable[[int, int], bytes], *, byte_budget: int = _MAX_DIRECTORY_BYTES
) -> tuple[Mapping[str, float | int | str], int]:
    """Parse a TIFF via small targeted ranges, returning facts and bytes requested.

    TIFF permits its first image-file directory to occur far from byte zero.
    Fetching that directory and the two rational scale values by ``Range`` keeps
    the inspection bounded even for an enormous whole-slide image.
    """
    if byte_budget < 16:
        raise TiffHeaderError("TIFF byte budget must cover an initial BigTIFF header")
    requested = 0

    def read(offset: int, length: int) -> bytes:
        nonlocal requested
        if offset < 0 or length < 0 or requested + length > byte_budget:
            raise TiffHeaderError("needed TIFF header regions exceed the fixed byte budget")
        raw = read_range(offset, length)
        if len(raw) != length:
            raise TiffHeaderError("TIFF range response was shorter than requested")
        requested += length
        return raw

    initial = read(0, 16)
    (
        order,
        ifd_offset,
        count_size,
        count_field_size,
        entry_size,
        inline_size,
        count_format,
        offset_format,
        kind,
    ) = _header_layout(initial)
    entry_count = struct.unpack(f"{order}{count_format}", read(ifd_offset, count_size))[0]
    directory_size = entry_count * entry_size
    if directory_size > _MAX_DIRECTORY_BYTES:
        raise TiffHeaderError("first TIFF directory exceeds the fixed directory-byte limit")
    directory = read(ifd_offset + count_size, directory_size)

    tags: dict[int, int | tuple[int, int]] = {}
    for index in range(entry_count):
        position = index * entry_size
        tag, field_type = struct.unpack_from(f"{order}HH", directory, position)
        count = struct.unpack_from(f"{order}{count_format}", directory, position + 4)[0]
        value_field = position + 4 + count_field_size
        value_offset = struct.unpack_from(f"{order}{offset_format}", directory, value_field)[0]
        if tag in _REQUIRED_TAGS:
            size = _TYPE_SIZES.get(field_type)
            if size is None or count != 1:
                raise TiffHeaderError(
                    f"unsupported TIFF tag type/count: type={field_type}, count={count}"
                )
            raw = (
                directory[value_field : value_field + size]
                if size <= inline_size
                else read(value_offset, size)
            )
            value = _decode_value(raw, field_type=field_type, count=count, order=order)
            tags[tag] = value if field_type == 5 else value[0]

    missing = sorted(_REQUIRED_TAGS - tags.keys())
    if missing:
        raise TiffHeaderError(f"TIFF header lacks required tags {missing}")
    x_num, x_den = tags[_X_RESOLUTION]  # type: ignore[misc]
    y_num, y_den = tags[_Y_RESOLUTION]  # type: ignore[misc]
    if not all((x_num, x_den, y_num, y_den)):
        raise TiffHeaderError("TIFF header has a zero physical-resolution value")
    unit = int(tags[_RESOLUTION_UNIT])
    microns_per_unit = {2: 25_400.0, 3: 10_000.0}.get(unit)
    if microns_per_unit is None:
        raise TiffHeaderError(f"TIFF resolution unit {unit} cannot be converted to microns")
    return (
        {
            "tiff_kind": kind,
            "pixel_width": int(tags[_PIXEL_WIDTH]),
            "pixel_height": int(tags[_PIXEL_HEIGHT]),
            "resolution_unit": {2: "inch", 3: "centimeter"}[unit],
            "x_microns_per_pixel": microns_per_unit * x_den / x_num,
            "y_microns_per_pixel": microns_per_unit * y_den / y_num,
        },
        requested,
    )


def parse_tiff_header(data: bytes) -> Mapping[str, float | int | str]:
    """Extract first-page pixel dimensions and physical sampling from TIFF bytes."""
    if len(data) < 8:
        raise TiffHeaderError("header prefix is shorter than a TIFF header")
    if data[:2] not in {b"II", b"MM"}:
        raise TiffHeaderError("not a TIFF byte-order marker")

    def read_range(offset: int, length: int) -> bytes:
        return data[offset : offset + length]

    # This in-memory compatibility wrapper can revisit a few ranges. The remote
    # path, by contrast, has a fixed 64 KiB cumulative budget.
    result, _ = parse_tiff_header_ranges(read_range, byte_budget=max(len(data) * 3, 16))
    return result
