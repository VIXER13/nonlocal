from __future__ import annotations

import pytest

from textual.app import App, ComposeResult
from textual.widgets import Select

from nonloccfg.widgets.task_tab import TaskSection, TimeSection
from nonloccfg.widgets.options_tab import SaveSection
from nonloccfg.widgets.boundaries import ThermalBCWidget, MechanicalBCWidget


# ---------------------------------------------------------------------------
# Helpers — minimal host apps that mount a single widget for testing
# ---------------------------------------------------------------------------

def _app_for(widget_factory):
    """Return a one-shot App class that composes a single widget."""
    class _Host(App):
        def compose(self) -> ComposeResult:
            yield widget_factory()
    return _Host


# ---------------------------------------------------------------------------
# TaskSection
# ---------------------------------------------------------------------------

class TestTaskSection:
    async def test_defaults(self):
        async with _app_for(TaskSection)().run_test() as pilot:
            section = pilot.app.query_one(TaskSection)
            data = section.collect()
            assert data["dimension"] == 2
            assert data["problem"] == "thermal"
            assert data["time_dependency"] is False

    async def test_load_roundtrip(self):
        async with _app_for(TaskSection)().run_test() as pilot:
            section = pilot.app.query_one(TaskSection)
            original = {"problem": "mechanical", "time_dependency": True}
            section.load(original)
            await pilot.pause()
            result = section.collect()
            assert result["problem"] == "mechanical"
            assert result["time_dependency"] is True

    async def test_thermomechanical_roundtrip(self):
        async with _app_for(TaskSection)().run_test() as pilot:
            section = pilot.app.query_one(TaskSection)
            section.load({"problem": "thermomechanical", "time_dependency": False})
            await pilot.pause()
            assert section.collect()["problem"] == "thermomechanical"


# ---------------------------------------------------------------------------
# TimeSection
# ---------------------------------------------------------------------------

class TestTimeSection:
    async def test_empty_collect(self):
        async with _app_for(TimeSection)().run_test() as pilot:
            data = pilot.app.query_one(TimeSection).collect()
            assert data == {}

    async def test_load_roundtrip(self):
        async with _app_for(TimeSection)().run_test() as pilot:
            section = pilot.app.query_one(TimeSection)
            original = {"time_step": 0.01, "steps_count": 100, "save_frequency": 5}
            section.load(original)
            await pilot.pause()
            result = section.collect()
            assert result["time_step"] == pytest.approx(0.01)
            assert result["steps_count"] == 100
            assert result["save_frequency"] == 5

    async def test_partial_load(self):
        async with _app_for(TimeSection)().run_test() as pilot:
            section = pilot.app.query_one(TimeSection)
            section.load({"time_step": 1.0})
            await pilot.pause()
            result = section.collect()
            assert result["time_step"] == pytest.approx(1.0)
            assert "steps_count" not in result


# ---------------------------------------------------------------------------
# SaveSection
# ---------------------------------------------------------------------------

class TestSaveSection:
    async def test_empty_collect(self):
        async with _app_for(SaveSection)().run_test() as pilot:
            assert pilot.app.query_one(SaveSection).collect() == {}

    async def test_load_roundtrip(self):
        async with _app_for(SaveSection)().run_test() as pilot:
            section = pilot.app.query_one(SaveSection)
            original = {"folder": "./out", "csv": "solution", "precision": 7}
            section.load(original)
            await pilot.pause()
            result = section.collect()
            assert result["folder"] == "./out"
            assert result["csv"] == "solution"
            assert result["precision"] == 7


# ---------------------------------------------------------------------------
# ThermalBCWidget
# ---------------------------------------------------------------------------

def _thermal_bc_app():
    class _Host(App):
        def compose(self) -> ComposeResult:
            yield ThermalBCWidget(0, "wall")
    return _Host()


class TestThermalBCWidget:
    async def test_default_kind_is_temperature(self):
        async with _thermal_bc_app().run_test() as pilot:
            data = pilot.app.query_one(ThermalBCWidget).collect()
            assert data["kind"] == "temperature"

    async def test_load_temperature_roundtrip(self):
        async with _thermal_bc_app().run_test() as pilot:
            w = pilot.app.query_one(ThermalBCWidget)
            w.load({"kind": "temperature", "temperature": 300.0})
            await pilot.pause()
            result = w.collect()
            assert result["kind"] == "temperature"
            assert result["temperature"] == pytest.approx(300.0)

    async def test_load_flux_roundtrip(self):
        async with _thermal_bc_app().run_test() as pilot:
            w = pilot.app.query_one(ThermalBCWidget)
            w.load({"kind": "flux", "flux": 500.0})
            await pilot.pause()
            result = w.collect()
            assert result["kind"] == "flux"
            assert result["flux"] == pytest.approx(500.0)

    async def test_load_convection_roundtrip(self):
        async with _thermal_bc_app().run_test() as pilot:
            w = pilot.app.query_one(ThermalBCWidget)
            w.load({"kind": "convection", "temperature": 293.0, "heat_transfer_coefficient": 25.0})
            await pilot.pause()
            result = w.collect()
            assert result["kind"] == "convection"
            assert result["temperature"] == pytest.approx(293.0)
            assert result["heat_transfer_coefficient"] == pytest.approx(25.0)


# ---------------------------------------------------------------------------
# MechanicalBCWidget
# ---------------------------------------------------------------------------

def _mech_bc_app():
    class _Host(App):
        def compose(self) -> ComposeResult:
            yield MechanicalBCWidget(0, "support")
    return _Host()


class TestMechanicalBCWidget:
    async def test_default_collect_is_null_null(self):
        async with _mech_bc_app().run_test() as pilot:
            result = pilot.app.query_one(MechanicalBCWidget).collect()
            assert result == [None, None]

    async def test_load_displacement_x(self):
        async with _mech_bc_app().run_test() as pilot:
            w = pilot.app.query_one(MechanicalBCWidget)
            w.load([{"displacement": 0.0}, None])
            await pilot.pause()
            result = w.collect()
            assert result[0] == {"displacement": 0.0}
            assert result[1] is None

    async def test_load_both_axes(self):
        async with _mech_bc_app().run_test() as pilot:
            w = pilot.app.query_one(MechanicalBCWidget)
            w.load([{"displacement": 0.0}, {"displacement": 0.0}])
            await pilot.pause()
            result = w.collect()
            assert result[0] == {"displacement": 0.0}
            assert result[1] == {"displacement": 0.0}
