from __future__ import annotations

import struct

import pytest

from src.reference.he_tiff_header import (
    TiffHeaderError,
    classify_resolution,
    parse_tiff_header,
    parse_tiff_header_ranges,
)


def _classic_tiff(*, include_unit: bool = True, resolution: int = 101_600) -> bytes:
    entries = [(256, 4, 1, 1000), (257, 4, 1, 500)]
    if include_unit:
        entries.append((296, 3, 1, 2))
    rational_start = 8 + 2 + (len(entries) + 2) * 12 + 4
    entries.extend([(282, 5, 1, rational_start), (283, 5, 1, rational_start + 8)])
    data = bytearray(b"II" + struct.pack("<H", 42) + struct.pack("<I", 8))
    data += struct.pack("<H", len(entries))
    for tag, field_type, count, value in entries:
        data += struct.pack("<HHI", tag, field_type, count)
        data += (
            struct.pack("<H", value) + b"\x00\x00" if field_type == 3 else struct.pack("<I", value)
        )
    data += b"\x00" * 4 + struct.pack("<II", resolution, 1) + struct.pack("<II", resolution, 1)
    return bytes(data)


def test_classic_tiff_header_reports_native_dimensions_and_scale():
    header = parse_tiff_header(_classic_tiff())
    assert header == {
        "tiff_kind": "classic_tiff",
        "pixel_width": 1000,
        "pixel_height": 500,
        "resolution_unit": "inch",
        "x_resolution_numerator": 101_600,
        "x_resolution_denominator": 1,
        "y_resolution_numerator": 101_600,
        "y_resolution_denominator": 1,
        "resolution_status": "measured",
        "resolution_not_measured_reason": "",
        "x_microns_per_pixel": pytest.approx(0.25),
        "y_microns_per_pixel": pytest.approx(0.25),
    }


def test_a_writer_default_resolution_is_not_reported_as_a_measured_scale():
    """72 dpi is what a writer emits when it knows nothing about the specimen."""
    header = parse_tiff_header(_classic_tiff(resolution=72))
    assert header["resolution_status"] == "writer_default"
    assert header["x_microns_per_pixel"] is None
    assert header["y_microns_per_pixel"] is None
    # The raw tag survives, so the classification can be checked from the table.
    assert header["x_resolution_numerator"] == 72
    assert "writer default" in header["resolution_not_measured_reason"]


def test_a_resolution_no_scanner_could_produce_is_refused_as_a_scale():
    header = parse_tiff_header(_classic_tiff(resolution=300))
    assert header["resolution_status"] == "implausible"
    assert header["x_microns_per_pixel"] is None
    assert "plausibility bound" in header["resolution_not_measured_reason"]


def test_the_htan_observed_default_resolves_to_no_scale():
    """The exact value every parsed HTA11 H&E header carried: 25400/72."""
    status, reason = classify_resolution(72.0, 2, 25_400.0 / 72.0)
    assert status == "writer_default"
    assert reason


def test_required_scale_tag_missing_from_header_is_a_refusal_not_a_default():
    with pytest.raises(TiffHeaderError, match="lacks required tags"):
        parse_tiff_header(_classic_tiff(include_unit=False))


def test_not_a_tiff_is_never_treated_as_a_native_image_header():
    with pytest.raises(TiffHeaderError, match="byte-order marker"):
        parse_tiff_header(b"not-a-tiff")


def test_targeted_ranges_parse_a_directory_that_is_not_in_the_initial_prefix():
    ifd_offset = 100_000
    entries = [(256, 4, 1, 1000), (257, 4, 1, 500), (296, 3, 1, 2)]
    rational_start = ifd_offset + 2 + (len(entries) + 2) * 12 + 4
    entries.extend([(282, 5, 1, rational_start), (283, 5, 1, rational_start + 8)])
    image = bytearray(b"II" + struct.pack("<H", 42) + struct.pack("<I", ifd_offset))
    image.extend(b"\x00" * (ifd_offset - len(image)))
    image += struct.pack("<H", len(entries))
    for tag, field_type, count, value in entries:
        image += struct.pack("<HHI", tag, field_type, count)
        inline_value = (
            struct.pack("<H", value) + b"\x00\x00"
            if field_type == 3
            else struct.pack("<I", value)
        )
        image += inline_value
    image += b"\x00" * 4 + struct.pack("<II", 101_600, 1) + struct.pack("<II", 101_600, 1)
    calls: list[tuple[int, int]] = []

    def read_range(offset: int, length: int) -> bytes:
        calls.append((offset, length))
        return bytes(image[offset : offset + length])

    facts, bytes_requested = parse_tiff_header_ranges(read_range, byte_budget=512)
    assert facts["pixel_width"] == 1000
    assert facts["x_microns_per_pixel"] == pytest.approx(0.25)
    assert bytes_requested < 512
    assert any(offset == ifd_offset for offset, _ in calls)
