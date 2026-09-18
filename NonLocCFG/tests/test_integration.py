"""
Integration tests — drive NonlocCfgApp through Textual's run_test() Pilot.

Each test boots the full app (all three tabs) and exercises real user
interactions: tab switching, problem-type changes, config dump/load,
Materials DB edits.  No mocking of internal widgets.
"""
from __future__ import annotations

import json
import textwrap

import pytest

from nonloccfg.app import NonlocCfgApp
from nonloccfg.materials_db import MaterialsDB
from nonloccfg.widgets.materials_db_tab import MaterialsDBTab
from nonloccfg.widgets.options_tab import OptionsTab
from nonloccfg.widgets.task_tab import TaskSection, TaskTab, TimeSection
from textual.widgets import Input, Select, TabbedContent


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def su2_mesh(tmp_path):
    """Minimal SU2 file with one boundary and one material marker."""
    content = textwrap.dedent("""\
        MARKER_TAG=Wall
        MARKER_ELEMS=1
        3 0 1
        MARKER_TAG=Body
        MARKER_ELEMS=1
        5 0 1 2
    """)
    path = tmp_path / "mesh.su2"
    path.write_text(content, encoding="utf-8")
    return path


@pytest.fixture
def thermal_config(tmp_path, su2_mesh):
    """Minimal thermal JSON config pointing at the temp mesh."""
    cfg = {
        "task": {"dimension": 2, "problem": "thermal", "time_dependency": False},
        "mesh": {"path": str(su2_mesh)},
        "boundaries": {"Wall": {"kind": "temperature", "temperature": 300.0}},
        "materials": {"Body": {"physical": {"conductivity": 50.0, "density": 7800.0}}},
    }
    path = tmp_path / "thermal.json"
    path.write_text(json.dumps(cfg, indent=2), encoding="utf-8")
    return path


@pytest.fixture
def thermo_config(tmp_path, su2_mesh):
    """Thermomechanical config."""
    cfg = {
        "task": {"dimension": 2, "problem": "thermomechanical", "time_dependency": False},
        "mesh": {"path": str(su2_mesh)},
        "thermal_boundaries": {"Wall": {"kind": "temperature", "temperature": 300.0}},
        "mechanical_boundaries": {},
        "materials": {
            "Body": {
                "physical": {
                    "conductivity": 50.0,
                    "density": 7800.0,
                    "youngs_modulus": 200e9,
                    "poissons_ratio": 0.3,
                },
                "thermal_model": {"local_weight": 0.5, "nonlocal_radius": 0.1},
                "mechanical_model": {"local_weight": 0.5, "nonlocal_radius": 0.1},
            }
        },
    }
    path = tmp_path / "thermo.json"
    path.write_text(json.dumps(cfg, indent=2), encoding="utf-8")
    return path


@pytest.fixture
def isolated_db(tmp_path, monkeypatch):
    """Replace the module-level `db` singleton with an empty temp-dir DB."""
    import nonloccfg.db as db_module
    import nonloccfg.widgets.materials_db_tab as tab_module
    import nonloccfg.widgets.material_picker as picker_module

    fresh = MaterialsDB(tmp_path / "test_db.json")
    monkeypatch.setattr(db_module, "db", fresh)
    monkeypatch.setattr(tab_module, "db", fresh)
    monkeypatch.setattr(picker_module, "db", fresh)
    return fresh


# ---------------------------------------------------------------------------
# App startup
# ---------------------------------------------------------------------------

class TestAppStartup:
    async def test_launches_without_error(self):
        async with NonlocCfgApp().run_test() as pilot:
            await pilot.pause()
            assert pilot.app.query_one(TaskTab) is not None

    async def test_all_three_tabs_present(self):
        async with NonlocCfgApp().run_test() as pilot:
            tc = pilot.app.query_one(TabbedContent)
            tab_ids = {pane.id for pane in tc.query("TabPane")}
            assert {"tab-task", "tab-materials", "tab-options"} <= tab_ids

    async def test_default_problem_is_thermal(self):
        async with NonlocCfgApp().run_test() as pilot:
            sel = pilot.app.query_one("#task-problem", Select)
            assert str(sel.value) == "thermal"


# ---------------------------------------------------------------------------
# Tab navigation
# ---------------------------------------------------------------------------

class TestTabNavigation:
    async def test_switch_to_materials_tab(self):
        async with NonlocCfgApp().run_test() as pilot:
            await pilot.click("#--content-tab-tab-materials")
            await pilot.pause()
            tc = pilot.app.query_one(TabbedContent)
            assert tc.active == "tab-materials"

    async def test_switch_to_options_tab(self):
        async with NonlocCfgApp().run_test() as pilot:
            await pilot.click("#--content-tab-tab-options")
            await pilot.pause()
            assert pilot.app.query_one(TabbedContent).active == "tab-options"

    async def test_switch_back_to_task_tab(self):
        async with NonlocCfgApp().run_test() as pilot:
            await pilot.click("#--content-tab-tab-options")
            await pilot.pause()
            await pilot.click("#--content-tab-tab-task")
            await pilot.pause()
            assert pilot.app.query_one(TabbedContent).active == "tab-task"


# ---------------------------------------------------------------------------
# Task tab — problem type changes
# ---------------------------------------------------------------------------

class TestProblemTypeChange:
    async def test_time_section_hidden_by_default(self):
        async with NonlocCfgApp().run_test() as pilot:
            assert not pilot.app.query_one(TimeSection).display

    async def test_time_section_shows_when_time_dependency_true(self):
        async with NonlocCfgApp().run_test() as pilot:
            # press "true" radio (second button in task-time-dependency)
            await pilot.click("#task-time-dependency RadioButton:last-of-type")
            await pilot.pause()
            assert pilot.app.query_one(TimeSection).display

    async def test_change_problem_to_mechanical(self):
        async with NonlocCfgApp().run_test() as pilot:
            sel = pilot.app.query_one("#task-problem", Select)
            sel.value = "mechanical"
            await pilot.pause()
            cfg = pilot.app.query_one(TaskTab).build_config()
            assert cfg["task"]["problem"] == "mechanical"

    async def test_change_problem_to_thermomechanical(self):
        async with NonlocCfgApp().run_test() as pilot:
            sel = pilot.app.query_one("#task-problem", Select)
            sel.value = "thermomechanical"
            await pilot.pause()
            cfg = pilot.app.query_one(TaskTab).build_config()
            assert cfg["task"]["problem"] == "thermomechanical"


# ---------------------------------------------------------------------------
# Config dump (Ctrl+S)
# ---------------------------------------------------------------------------

class TestConfigDump:
    async def test_dump_creates_file(self, tmp_path):
        out = tmp_path / "out.json"
        async with NonlocCfgApp().run_test() as pilot:
            pilot.app.query_one(OptionsTab).query_one("#output-path", Input).value = str(out)
            await pilot.press("ctrl+s")
            await pilot.pause()
        assert out.exists()

    async def test_dump_is_valid_json(self, tmp_path):
        out = tmp_path / "out.json"
        async with NonlocCfgApp().run_test() as pilot:
            pilot.app.query_one(OptionsTab).query_one("#output-path", Input).value = str(out)
            await pilot.press("ctrl+s")
            await pilot.pause()
        cfg = json.loads(out.read_text())
        assert "task" in cfg
        assert cfg["task"]["problem"] == "thermal"

    async def test_dump_includes_save_section(self, tmp_path):
        out = tmp_path / "out.json"
        async with NonlocCfgApp().run_test() as pilot:
            # fill in save.folder
            pilot.app.query_one("#save-folder", Input).value = "./results"
            pilot.app.query_one(OptionsTab).query_one("#output-path", Input).value = str(out)
            await pilot.press("ctrl+s")
            await pilot.pause()
        cfg = json.loads(out.read_text())
        assert cfg.get("save", {}).get("folder") == "./results"


# ---------------------------------------------------------------------------
# Config load round-trip
# ---------------------------------------------------------------------------

class TestConfigLoadRoundtrip:
    async def test_thermal_config_loads(self, thermal_config):
        async with NonlocCfgApp().run_test() as pilot:
            app = pilot.app
            cfg = json.loads(thermal_config.read_text())
            app.query_one(OptionsTab).load_options(cfg)
            await app.run_worker(
                app.query_one(TaskTab).load_config(cfg, thermal_config.parent),
                exclusive=True,
            ).wait()
            await pilot.pause()
            result = app.query_one(TaskTab).build_config()
            assert result["task"]["problem"] == "thermal"
            assert result["task"]["time_dependency"] is False

    async def test_thermal_config_boundaries_loaded(self, thermal_config):
        async with NonlocCfgApp().run_test() as pilot:
            app = pilot.app
            cfg = json.loads(thermal_config.read_text())
            await app.run_worker(
                app.query_one(TaskTab).load_config(cfg, thermal_config.parent),
                exclusive=True,
            ).wait()
            await pilot.pause()
            result = app.query_one(TaskTab).build_config()
            assert "Wall" in result.get("boundaries", {})

    async def test_thermal_config_materials_loaded(self, thermal_config):
        async with NonlocCfgApp().run_test() as pilot:
            app = pilot.app
            cfg = json.loads(thermal_config.read_text())
            await app.run_worker(
                app.query_one(TaskTab).load_config(cfg, thermal_config.parent),
                exclusive=True,
            ).wait()
            await pilot.pause()
            result = app.query_one(TaskTab).build_config()
            assert "Body" in result.get("materials", {})

    async def test_thermomechanical_config_loads(self, thermo_config):
        async with NonlocCfgApp().run_test() as pilot:
            app = pilot.app
            cfg = json.loads(thermo_config.read_text())
            await app.run_worker(
                app.query_one(TaskTab).load_config(cfg, thermo_config.parent),
                exclusive=True,
            ).wait()
            await pilot.pause()
            result = app.query_one(TaskTab).build_config()
            assert result["task"]["problem"] == "thermomechanical"
            mat = result["materials"]["Body"]
            assert "thermal_model" in mat
            assert "mechanical_model" in mat

    async def test_dump_load_dump_produces_same_task(self, tmp_path, thermal_config):
        """Load a config, dump it, reload the dump — task section must match."""
        out = tmp_path / "roundtrip.json"
        async with NonlocCfgApp().run_test() as pilot:
            app = pilot.app
            cfg = json.loads(thermal_config.read_text())
            await app.run_worker(
                app.query_one(TaskTab).load_config(cfg, thermal_config.parent),
                exclusive=True,
            ).wait()
            app.query_one(OptionsTab).query_one("#output-path", Input).value = str(out)
            await pilot.press("ctrl+s")
            await pilot.pause()

        dumped = json.loads(out.read_text())
        assert dumped["task"] == cfg["task"]


# ---------------------------------------------------------------------------
# Materials DB tab
# ---------------------------------------------------------------------------

class TestMaterialsDBTab:
    async def test_new_material_saves(self, isolated_db):
        # Use a tall terminal so #btn-save (inside VerticalScroll) is visible
        async with NonlocCfgApp().run_test(size=(160, 60)) as pilot:
            await pilot.click("#--content-tab-tab-materials")
            await pilot.pause()
            await pilot.click("#btn-new")
            await pilot.pause()
            inp = pilot.app.query_one("#field-name", Input)
            inp.value = "TestMat"
            await pilot.pause()
            await pilot.click("#btn-save")
            await pilot.pause()
        assert "TestMat" in isolated_db.all_names()

    async def test_delete_material(self, isolated_db):
        isolated_db.upsert("ToDelete", {"density": 1000.0})
        async with NonlocCfgApp().run_test(size=(160, 60)) as pilot:
            await pilot.click("#--content-tab-tab-materials")
            await pilot.pause()
            lv = pilot.app.query_one("#db-list")
            lv.index = 0
            await pilot.pause()
            await pilot.click("#btn-delete")
            await pilot.pause()
        assert "ToDelete" not in isolated_db.all_names()

    async def test_edit_existing_material(self, isolated_db):
        isolated_db.upsert("Iron", {"density": 7800.0})
        async with NonlocCfgApp().run_test(size=(160, 60)) as pilot:
            await pilot.click("#--content-tab-tab-materials")
            await pilot.pause()
            lv = pilot.app.query_one("#db-list")
            lv.index = 0
            await pilot.pause()
            tab = pilot.app.query_one(MaterialsDBTab)
            tab.query_one("#field-density", Input).value = "7874.0"
            await pilot.pause()
            await pilot.click("#btn-save")
            await pilot.pause()
        assert isolated_db.get("Iron")["density"] == pytest.approx(7874.0)

    async def test_save_without_name_shows_no_entry(self, isolated_db):
        async with NonlocCfgApp().run_test(size=(160, 60)) as pilot:
            await pilot.click("#--content-tab-tab-materials")
            await pilot.pause()
            await pilot.click("#btn-new")
            await pilot.pause()
            # Leave name blank, try to save
            await pilot.click("#btn-save")
            await pilot.pause()
        assert isolated_db.all_names() == []
