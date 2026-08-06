from __future__ import annotations

import json
import pytest

from nonloccfg.materials_db import MaterialsDB, _fuzzy_score


# ---------------------------------------------------------------------------
# _fuzzy_score
# ---------------------------------------------------------------------------

class TestFuzzyScore:
    def test_empty_query_always_matches(self):
        assert _fuzzy_score("", "Steel") > 0

    def test_exact_substring_scores_high(self):
        assert _fuzzy_score("steel", "Steel") >= 50

    def test_substring_at_start_scores_higher_than_end(self):
        start = _fuzzy_score("al", "Aluminum")
        end   = _fuzzy_score("um", "Aluminum")
        assert start > end

    def test_no_match_returns_zero(self):
        assert _fuzzy_score("zzz", "Steel") == 0

    def test_subsequence_matches(self):
        # "stl" is a subsequence of "Steel"
        assert _fuzzy_score("stl", "Steel") > 0

    def test_case_insensitive(self):
        assert _fuzzy_score("STEEL", "steel") > 0


# ---------------------------------------------------------------------------
# MaterialsDB
# ---------------------------------------------------------------------------

@pytest.fixture
def db(tmp_path):
    return MaterialsDB(tmp_path / "test_db.json")


class TestMaterialsDB:
    def test_empty_db_has_no_names(self, db):
        assert db.all_names() == []

    def test_upsert_and_get(self, db):
        db.upsert("Steel", {"density": 7800.0})
        assert db.get("Steel") == {"density": 7800.0}

    def test_get_missing_returns_none(self, db):
        assert db.get("Nonexistent") is None

    def test_all_names_sorted(self, db):
        db.upsert("Zinc", {})
        db.upsert("Aluminum", {})
        db.upsert("Copper", {})
        assert db.all_names() == ["Aluminum", "Copper", "Zinc"]

    def test_upsert_filters_none_values(self, db):
        db.upsert("Mat", {"density": 1000.0, "capacity": None})
        assert "capacity" not in db.get("Mat")

    def test_upsert_overwrites_existing(self, db):
        db.upsert("Steel", {"density": 7800.0})
        db.upsert("Steel", {"density": 7850.0})
        assert db.get("Steel")["density"] == 7850.0

    def test_delete_removes_entry(self, db):
        db.upsert("Steel", {"density": 7800.0})
        db.delete("Steel")
        assert db.get("Steel") is None
        assert "Steel" not in db.all_names()

    def test_delete_nonexistent_is_noop(self, db):
        db.delete("Ghost")  # must not raise

    def test_flush_persists_to_disk(self, tmp_path):
        path = tmp_path / "db.json"
        db1 = MaterialsDB(path)
        db1.upsert("Iron", {"density": 7874.0})

        db2 = MaterialsDB(path)
        assert db2.get("Iron") == {"density": 7874.0}

    def test_corrupted_file_loads_empty(self, tmp_path):
        path = tmp_path / "db.json"
        path.write_text("not json", encoding="utf-8")
        db = MaterialsDB(path)
        assert db.all_names() == []

    def test_get_returns_copy(self, db):
        db.upsert("Steel", {"density": 7800.0})
        copy = db.get("Steel")
        copy["density"] = 0
        assert db.get("Steel")["density"] == 7800.0

    def test_search_returns_matches(self, db):
        db.upsert("Steel", {})
        db.upsert("Aluminum", {})
        db.upsert("Copper", {})
        results = db.search("al")
        assert "Aluminum" in results

    def test_search_excludes_non_matches(self, db):
        db.upsert("Steel", {})
        db.upsert("Copper", {})
        results = db.search("zzz")
        assert results == []

    def test_search_orders_by_score(self, db):
        db.upsert("Aluminum", {})
        db.upsert("Alloy Steel", {})
        results = db.search("al")
        # Both contain "al"; the one where it starts earlier scores higher
        assert results[0] == "Aluminum"
