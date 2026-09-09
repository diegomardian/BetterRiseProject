"""Read native TIFF header facts without downloading an H&E image.

Only the beginning of a TIFF is needed to find its first image-file directory
and its pixel-size tags in the usual baseline-TIFF layout.  This module makes
no claim that a deposited image is the scanner's original full-resolution
object; it reports only what its native TIFF header declares.
"""

from __future__ import annotations

import struct
from collections.abc import Mapping


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


def _read_value(
    data: bytes,
    *,
    offset: int,
    field_type: int,
    count: int,
    value_offset: int,
    inline_size: int,
    order: str,
) -> int | tuple[int, int]:
    if count != 1 or field_type not in _TYPE_SIZES:
        raise TiffHeaderError(f"unsupported TIFF tag type/count: type={field_type}, count={count}")
    size = _TYPE_SIZES[field_type]
    start = value_offset if size > inline_size else offset
    end = start + size
    if start < 0 or end > len(data):
        raise TiffHeaderError("needed TIFF tag is outside the bounded header prefix")
    return struct.unpack_from(f"{order}{_TYPE_FORMATS[field_type]}", data, start)


def parse_tiff_header(data: bytes) -> Mapping[str, float | int | str]:
    """Extract first-page pixel dimensions and physical sampling from TIFF bytes."""
    if len(data) < 8:
        raise TiffHeaderError("header prefix is shorter than a TIFF header")
    endian = data[:2]
    if endian == b"II":
        order = "<"
    elif endian == b"MM":
        order = ">"
    else:
        raise TiffHeaderError("not a TIFF byte-order marker")
    magic = struct.unpack_from(f"{order}H", data, 2)[0]
    if magic == 42:
        ifd_offset = struct.unpack_from(f"{order}I", data, 4)[0]
        count_size, count_field_size, entry_size, inline_size = 2, 4, 12, 4
        count_format, offset_format = "H", "I"
        kind = "classic_tiff"
    elif magic == 43:
        if len(data) < 16:
            raise TiffHeaderError("header prefix is shorter than a BigTIFF header")
        offset_size, reserved = struct.unpack_from(f"{order}HH", data, 4)
        if offset_size != 8 or reserved != 0:
            raise TiffHeaderError("unsupported BigTIFF offset layout")
        ifd_offset = struct.unpack_from(f"{order}Q", data, 8)[0]
        count_size, count_field_size, entry_size, inline_size = 8, 8, 20, 8
        count_format, offset_format = "Q", "Q"
        kind = "big_tiff"
    else:
        raise TiffHeaderError(f"unsupported TIFF magic value {magic}")

    if ifd_offset + count_size > len(data):
        raise TiffHeaderError("first TIFF directory is outside the bounded header prefix")
    entry_count = struct.unpack_from(f"{order}{count_format}", data, ifd_offset)[0]
    entries_start = ifd_offset + count_size
    if entries_start + entry_count * entry_size > len(data):
        raise TiffHeaderError("first TIFF directory is outside the bounded header prefix")

    tags: dict[int, int | tuple[int, int]] = {}
    for index in range(entry_count):
        position = entries_start + index * entry_size
        tag, field_type = struct.unpack_from(f"{order}HH", data, position)
        count = struct.unpack_from(f"{order}{count_format}", data, position + 4)[0]
        value_field = position + 4 + count_field_size
        value_offset = struct.unpack_from(f"{order}{offset_format}", data, value_field)[0]
        if tag in {_PIXEL_WIDTH, _PIXEL_HEIGHT, _X_RESOLUTION, _Y_RESOLUTION, _RESOLUTION_UNIT}:
            value = _read_value(
                data,
                offset=value_field,
                field_type=field_type,
                count=count,
                value_offset=value_offset,
                inline_size=inline_size,
                order=order,
            )
            tags[tag] = value if field_type == 5 else value[0]

    missing = sorted(
        {_PIXEL_WIDTH, _PIXEL_HEIGHT, _X_RESOLUTION, _Y_RESOLUTION, _RESOLUTION_UNIT} - tags.keys()
    )
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
    return {
        "tiff_kind": kind,
        "pixel_width": int(tags[_PIXEL_WIDTH]),
        "pixel_height": int(tags[_PIXEL_HEIGHT]),
        "resolution_unit": {2: "inch", 3: "centimeter"}[unit],
        "x_microns_per_pixel": microns_per_unit * x_den / x_num,
        "y_microns_per_pixel": microns_per_unit * y_den / y_num,
    }
