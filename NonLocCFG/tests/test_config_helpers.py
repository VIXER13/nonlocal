from __future__ import annotations

import pytest

from nonloccfg.config_helpers import _format_value, _parse_value


class TestParseValue:
    def test_empty_string_returns_none(self):
        assert _parse_value("") is None

    def test_whitespace_only_returns_none(self):
        assert _parse_value("   ") is None

    def test_integer(self):
        assert _parse_value("42") == 42
        assert isinstance(_parse_value("42"), int)

    def test_negative_integer(self):
        assert _parse_value("-5") == -5

    def test_float(self):
        result = _parse_value("3.14")
        assert isinstance(result, float)
        assert result == pytest.approx(3.14)

    def test_scientific_notation(self):
        result = _parse_value("2e9")
        assert isinstance(result, float)
        assert result == pytest.approx(2e9)

    def test_json_array_of_floats(self):
        assert _parse_value("[0.5, 0.5]") == [0.5, 0.5]

    def test_json_array_with_null(self):
        assert _parse_value("[null, 0.3]") == [None, 0.3]

    def test_plain_string(self):
        assert _parse_value("x y: expr") == "x y: expr"

    def test_whitespace_is_stripped(self):
        assert _parse_value("  10  ") == 10

    def test_invalid_json_array_falls_through_to_string(self):
        assert _parse_value("[not json") == "[not json"

    def test_zero_is_int(self):
        result = _parse_value("0")
        assert result == 0
        assert isinstance(result, int)


class TestFormatValue:
    def test_int(self):
        assert _format_value(42) == "42"

    def test_float(self):
        assert _format_value(3.14) == "3.14"

    def test_list_serialised_as_json(self):
        assert _format_value([0.5, 0.5]) == "[0.5, 0.5]"

    def test_list_with_null(self):
        assert _format_value([None, 0.3]) == "[null, 0.3]"

    def test_string_passthrough(self):
        assert _format_value("hello") == "hello"

    def test_roundtrip_float(self):
        assert _parse_value(_format_value(1.23)) == pytest.approx(1.23)

    def test_roundtrip_list(self):
        original = [1.0, 2.0]
        assert _parse_value(_format_value(original)) == original
