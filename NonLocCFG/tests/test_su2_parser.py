from __future__ import annotations

import textwrap

import pytest

from nonloccfg.su2_parser import parse_su2_markers


def _write_su2(tmp_path, content: str) -> str:
    path = tmp_path / "mesh.su2"
    path.write_text(textwrap.dedent(content), encoding="utf-8")
    return str(path)


class TestParseSu2Markers:
    def test_missing_file_returns_empty(self, tmp_path):
        btags, mtags = parse_su2_markers(str(tmp_path / "missing.su2"))
        assert btags == []
        assert mtags == []

    def test_unknown_dimension_returns_empty(self, tmp_path):
        path = _write_su2(tmp_path, "MARKER_TAG=wall\nMARKER_ELEMS=1\n3 0 1\n")
        btags, mtags = parse_su2_markers(path, dimension=99)
        assert btags == []
        assert mtags == []

    def test_2d_boundary_line_element(self, tmp_path):
        # elem type 3 = line → boundary in 2D
        content = """
            MARKER_TAG=wall
            MARKER_ELEMS=1
            3 0 1
        """
        path = _write_su2(tmp_path, content)
        btags, mtags = parse_su2_markers(path, dimension=2)
        assert btags == ["wall"]
        assert mtags == []

    def test_2d_material_triangle_element(self, tmp_path):
        # elem type 5 = triangle → material in 2D (not a boundary type)
        content = """
            MARKER_TAG=body
            MARKER_ELEMS=1
            5 0 1 2
        """
        path = _write_su2(tmp_path, content)
        btags, mtags = parse_su2_markers(path, dimension=2)
        assert btags == []
        assert mtags == ["body"]

    def test_multiple_markers(self, tmp_path):
        content = """
            MARKER_TAG=left_wall
            MARKER_ELEMS=2
            3 0 1
            3 1 2
            MARKER_TAG=right_wall
            MARKER_ELEMS=1
            3 3 4
            MARKER_TAG=fluid
            MARKER_ELEMS=1
            5 0 1 2
        """
        path = _write_su2(tmp_path, content)
        btags, mtags = parse_su2_markers(path, dimension=2)
        assert set(btags) == {"left_wall", "right_wall"}
        assert mtags == ["fluid"]

    def test_marker_classified_only_once(self, tmp_path):
        # Two element lines under same marker — tag should appear only once
        content = """
            MARKER_TAG=wall
            MARKER_ELEMS=2
            3 0 1
            3 2 3
        """
        path = _write_su2(tmp_path, content)
        btags, _ = parse_su2_markers(path, dimension=2)
        assert btags.count("wall") == 1

    def test_3d_boundary_triangle(self, tmp_path):
        # elem type 5 = triangle → boundary in 3D
        content = """
            MARKER_TAG=surface
            MARKER_ELEMS=1
            5 0 1 2
        """
        path = _write_su2(tmp_path, content)
        btags, mtags = parse_su2_markers(path, dimension=3)
        assert btags == ["surface"]
        assert mtags == []

    def test_1d_boundary_point(self, tmp_path):
        # elem type 15 = point → boundary in 1D
        content = """
            MARKER_TAG=endpoint
            MARKER_ELEMS=1
            15 0
        """
        path = _write_su2(tmp_path, content)
        btags, mtags = parse_su2_markers(path, dimension=1)
        assert btags == ["endpoint"]

    def test_empty_file_returns_empty(self, tmp_path):
        path = _write_su2(tmp_path, "")
        btags, mtags = parse_su2_markers(path)
        assert btags == []
        assert mtags == []
