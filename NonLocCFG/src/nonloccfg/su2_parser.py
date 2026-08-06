from __future__ import annotations

# SU2 element type codes by dimension of the element:
#   0D: 15=point
#   1D: 3=line
#   2D: 5=triangle, 9=quad
#   3D: 10=tet, 12=hex, 13=prism, 14=pyramid
_BOUNDARY_ELEM_TYPES: dict[int, frozenset[int]] = {
    1: frozenset({15}),   # point elements mark boundaries in 1D meshes
    2: frozenset({3}),    # line elements mark boundaries in 2D meshes
    3: frozenset({5, 9}), # tri/quad elements mark boundaries in 3D meshes
}


def parse_su2_markers(path: str, dimension: int = 2) -> tuple[list[str], list[str]]:
    """Return (boundary_tags, material_tags) from an SU2 file."""
    boundary_type_set = _BOUNDARY_ELEM_TYPES.get(dimension)
    if boundary_type_set is None:
        return [], []

    boundary_tags: list[str] = []
    material_tags: list[str] = []
    current_tag: str | None = None
    classified = False
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if line.startswith("MARKER_TAG="):
                    current_tag = line.split("=", 1)[1].strip()
                    classified = False
                elif line.startswith("MARKER_ELEMS="):
                    pass
                elif current_tag and not classified and line:
                    parts = line.split()
                    if parts and parts[0].isdigit():
                        elem_type = int(parts[0])
                        if elem_type in boundary_type_set:
                            boundary_tags.append(current_tag)
                        else:
                            material_tags.append(current_tag)
                        classified = True
    except OSError:
        pass
    return boundary_tags, material_tags
