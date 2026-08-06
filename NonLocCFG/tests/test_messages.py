from __future__ import annotations

from nonloccfg.messages import MeshLoaded, ProblemChanged, TagReleased


class TestProblemChanged:
    def test_stores_problem(self):
        msg = ProblemChanged("thermal")
        assert msg.problem == "thermal"


class TestMeshLoaded:
    def test_stores_tags(self):
        msg = MeshLoaded(["wall", "inlet"], ["body"])
        assert msg.boundary_tags == ["wall", "inlet"]
        assert msg.material_tags == ["body"]

    def test_empty_tags(self):
        msg = MeshLoaded([], [])
        assert msg.boundary_tags == []
        assert msg.material_tags == []


class TestTagReleased:
    def test_stores_tag_and_pool(self):
        msg = TagReleased("wall", "tbc")
        assert msg.tag == "wall"
        assert msg.pool == "tbc"
