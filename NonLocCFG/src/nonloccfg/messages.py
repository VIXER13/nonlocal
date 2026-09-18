from __future__ import annotations

from textual.message import Message


class ProblemChanged(Message):
    def __init__(self, problem: str) -> None:
        self.problem = problem
        super().__init__()


class MeshLoaded(Message):
    def __init__(self, boundary_tags: list[str], material_tags: list[str]) -> None:
        self.boundary_tags = boundary_tags
        self.material_tags = material_tags
        super().__init__()


class TagReleased(Message):
    """Posted by a BC/material widget when it is removed."""
    def __init__(self, tag: str, pool: str) -> None:
        self.tag = tag
        self.pool = pool  # "tbc" / "mbc" / "mat"
        super().__init__()
