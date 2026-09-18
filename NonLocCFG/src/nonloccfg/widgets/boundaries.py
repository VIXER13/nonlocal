from __future__ import annotations

from textual import on
from textual.app import ComposeResult
from textual.containers import Horizontal, Vertical
from textual.widget import Widget
from textual.widgets import Button, Input, Label, Select

from nonloccfg.config_helpers import _format_value, _input, _parse_value, _select
from nonloccfg.messages import TagReleased


class ThermalBCWidget(Widget):
    DEFAULT_CSS = """
    ThermalBCWidget { height: auto; padding-left: 2; margin-top: 1; }
    ThermalBCWidget .name-row { height: 1; align: left middle; }
    ThermalBCWidget .name-row Label { color: $text-muted; width: 2; }
    ThermalBCWidget .name-row Label.bc-tag { color: $text; width: auto; }
    ThermalBCWidget .name-row Button { width: 3; min-width: 3; }
    ThermalBCWidget .row { height: 1; align: left middle; padding-left: 2; margin-top: 1; }
    ThermalBCWidget .key { color: $text-muted; width: 18; }
    ThermalBCWidget .row Input, ThermalBCWidget .row Select { width: 1fr; }
    """

    def __init__(self, index: int, tag: str, prefix: str = "tbc") -> None:
        self._index = index
        self._tag = tag
        self._p = prefix
        super().__init__(id=f"{prefix}-{index}")

    def compose(self) -> ComposeResult:
        p, i = self._p, self._index
        with Horizontal(classes="name-row"):
            yield Label("- ")
            yield Label(self._tag, classes="bc-tag", id=f"{p}-tag-label-{i}")
            yield Button("×", variant="error", id=f"{p}-remove-{i}", compact=True)
        with Horizontal(classes="row"):
            yield Label("kind:", classes="key")
            yield Select(
                [("temperature", "temperature"), ("flux", "flux"), ("convection", "convection"),
                 ("radiation", "radiation"), ("combined", "combined")],
                value="temperature", id=f"{p}-kind-{i}", compact=True,
            )
        with Horizontal(classes="row", id=f"{p}-row-temperature-{i}"):
            yield Label("temperature:", classes="key")
            yield Input(placeholder="float or 'x y: expr'", id=f"{p}-temperature-{i}", compact=True)
        with Horizontal(classes="row", id=f"{p}-row-flux-{i}"):
            yield Label("flux:", classes="key")
            yield Input(placeholder="float or 'x y: expr'", id=f"{p}-flux-{i}", compact=True)
        with Horizontal(classes="row", id=f"{p}-row-heat-transfer-{i}"):
            yield Label("heat_transfer:", classes="key")
            yield Input(placeholder="float", id=f"{p}-heat-transfer-{i}", compact=True)
        with Horizontal(classes="row", id=f"{p}-row-emissivity-{i}"):
            yield Label("emissivity:", classes="key")
            yield Input(placeholder="float", id=f"{p}-emissivity-{i}", compact=True)

    def on_mount(self) -> None:
        p, i = self._p, self._index
        self.query_one(f"#{p}-row-flux-{i}").display = False
        self.query_one(f"#{p}-row-heat-transfer-{i}").display = False
        self.query_one(f"#{p}-row-emissivity-{i}").display = False

    def _apply_kind(self, kind: str) -> None:
        p, i = self._p, self._index
        self.query_one(f"#{p}-row-temperature-{i}").display = kind in ("temperature", "convection", "combined")
        self.query_one(f"#{p}-row-flux-{i}").display = kind in ("flux", "combined")
        self.query_one(f"#{p}-row-heat-transfer-{i}").display = kind in ("convection", "combined")
        self.query_one(f"#{p}-row-emissivity-{i}").display = kind in ("radiation", "combined")

    @on(Select.Changed)
    def _on_kind(self, event: Select.Changed) -> None:
        if event.select.id != f"{self._p}-kind-{self._index}":
            return
        self._apply_kind(str(event.value))

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == f"{self._p}-remove-{self._index}":
            event.stop()
            self.post_message(TagReleased(self._tag, self._p))
            self.remove()

    def collect(self) -> dict:
        p, i = self._p, self._index
        kind = _select(self, f"{p}-kind-{i}") or "temperature"
        d: dict = {"kind": kind}
        for field, wid in (
            ("temperature",              f"{p}-temperature-{i}"),
            ("flux",                     f"{p}-flux-{i}"),
            ("heat_transfer_coefficient", f"{p}-heat-transfer-{i}"),
            ("emissivity",               f"{p}-emissivity-{i}"),
        ):
            v = _parse_value(_input(self, wid))
            if v is not None:
                d[field] = v
        return d

    def load(self, data: dict) -> None:
        p, i = self._p, self._index
        kind = data.get("kind", "temperature")
        self.query_one(f"#{p}-kind-{i}", Select).value = kind
        self._apply_kind(kind)
        for field, wid in (
            ("temperature",              f"{p}-temperature-{i}"),
            ("flux",                     f"{p}-flux-{i}"),
            ("heat_transfer_coefficient", f"{p}-heat-transfer-{i}"),
            ("emissivity",               f"{p}-emissivity-{i}"),
        ):
            if field in data:
                self.query_one(f"#{wid}", Input).value = _format_value(data[field])


class MechanicalBCWidget(Widget):
    DEFAULT_CSS = """
    MechanicalBCWidget { height: auto; padding-left: 2; margin-top: 1; }
    MechanicalBCWidget .name-row { height: 1; align: left middle; }
    MechanicalBCWidget .name-row Label { color: $text-muted; width: 2; }
    MechanicalBCWidget .name-row Label.bc-tag { color: $text; width: auto; }
    MechanicalBCWidget .name-row Button { width: 3; min-width: 3; }
    MechanicalBCWidget .row { height: 1; align: left middle; padding-left: 2; margin-top: 1; }
    MechanicalBCWidget .key { color: $text-muted; width: 4; }
    MechanicalBCWidget .row Select { width: 24; }
    MechanicalBCWidget .row Input { width: 1fr; }
    """

    def __init__(self, index: int, tag: str, prefix: str = "mbc") -> None:
        self._index = index
        self._tag = tag
        self._p = prefix
        super().__init__(id=f"{prefix}-{index}")

    def compose(self) -> ComposeResult:
        p, i = self._p, self._index
        with Horizontal(classes="name-row"):
            yield Label("- ")
            yield Label(self._tag, classes="bc-tag")
            yield Button("×", variant="error", id=f"{p}-remove-{i}", compact=True)
        for axis in ("x", "y"):
            with Horizontal(classes="row"):
                yield Label(f"{axis}:", classes="key")
                yield Select(
                    [("displacement", "displacement"), ("pressure", "pressure")],
                    allow_blank=True, prompt="—  (null)",
                    id=f"{p}-type-{axis}-{i}", compact=True,
                )
                yield Input(
                    placeholder="float or 'x y: expr'",
                    id=f"{p}-val-{axis}-{i}", compact=True,
                )

    def on_mount(self) -> None:
        p, i = self._p, self._index
        for axis in ("x", "y"):
            self.query_one(f"#{p}-val-{axis}-{i}").display = False

    def _apply_axis(self, axis: str, kind: str) -> None:
        p, i = self._p, self._index
        self.query_one(f"#{p}-val-{axis}-{i}").display = bool(kind)

    @on(Select.Changed)
    def _on_axis(self, event: Select.Changed) -> None:
        p, i = self._p, self._index
        for axis in ("x", "y"):
            if event.select.id == f"{p}-type-{axis}-{i}":
                self._apply_axis(axis, "" if event.value is Select.NULL else str(event.value))
                break

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == f"{self._p}-remove-{self._index}":
            event.stop()
            self.post_message(TagReleased(self._tag, self._p))
            self.remove()

    def collect(self) -> list:
        p, i = self._p, self._index
        result = []
        for axis in ("x", "y"):
            kind = _select(self, f"{p}-type-{axis}-{i}")
            if not kind:
                result.append(None)
            else:
                v = _parse_value(_input(self, f"{p}-val-{axis}-{i}"))
                result.append({kind: v} if v is not None else {kind: 0})
        return result

    def load(self, data: list) -> None:
        p, i = self._p, self._index
        for axis, entry in zip(("x", "y"), data):
            if not entry:
                continue
            kind, val = next(iter(entry.items()))
            self.query_one(f"#{p}-type-{axis}-{i}", Select).value = kind
            self._apply_axis(axis, kind)
            if val is not None:
                self.query_one(f"#{p}-val-{axis}-{i}", Input).value = _format_value(val)


class BoundariesSection(Widget):
    DEFAULT_CSS = """
    BoundariesSection { height: auto; margin-bottom: 1; }
    BoundariesSection Vertical { height: auto; }
    BoundariesSection .add-row { height: 1; align: left middle; padding-left: 2; margin-top: 1; }
    BoundariesSection .add-row Label { color: $text-muted; width: 6; }
    BoundariesSection .add-row Select { width: 1fr; }
    """

    def __init__(self) -> None:
        super().__init__()
        self._tbc_available: list[str] = []
        self._mbc_available: list[str] = []
        self._next_tbc = 0
        self._next_mbc = 0
        self._problem = "thermal"

    def compose(self) -> ComposeResult:
        with Vertical(id="thermal-bc-group"):
            yield Label("[bold][#5f87ff]boundaries[/][/]:", id="thermal-bc-label")
            with Horizontal(classes="add-row", id="tbc-add-row"):
                yield Label("add:")
                yield Select([], allow_blank=True, prompt="— select boundary —",
                             id="tbc-add-select", compact=True)
        with Vertical(id="mechanical-bc-group"):
            yield Label("[bold][#5f87ff]mechanical_boundaries[/][/]:", id="mechanical-bc-label")
            with Horizontal(classes="add-row", id="mbc-add-row"):
                yield Label("add:")
                yield Select([], allow_blank=True, prompt="— select boundary —",
                             id="mbc-add-select", compact=True)

    def on_mount(self) -> None:
        self._apply_problem()
        self._disable_adds()

    def set_boundary_tags(self, tags: list[str]) -> None:
        self._tbc_available = list(tags)
        self._mbc_available = list(tags)
        self._refresh_selects()

    def _disable_adds(self) -> None:
        for sid in ("tbc-add-select", "mbc-add-select"):
            self.query_one(f"#{sid}", Select).disabled = True

    def _refresh_selects(self) -> None:
        for sid, pool in (("tbc-add-select", self._tbc_available),
                          ("mbc-add-select", self._mbc_available)):
            sel = self.query_one(f"#{sid}", Select)
            sel.set_options([(t, t) for t in pool])
            sel.disabled = len(pool) == 0

    @on(Select.Changed, "#tbc-add-select")
    def _add_tbc(self, event: Select.Changed) -> None:
        if event.value is Select.NULL:
            return
        tag = str(event.value)
        if tag not in self._tbc_available:
            return
        self._tbc_available.remove(tag)
        idx = self._next_tbc
        self._next_tbc += 1
        add_row = self.query_one("#tbc-add-row")
        self.query_one("#thermal-bc-group").mount(ThermalBCWidget(idx, tag, "tbc"), before=add_row)
        self.query_one("#tbc-add-select", Select).value = Select.NULL
        self._refresh_selects()

    @on(Select.Changed, "#mbc-add-select")
    def _add_mbc(self, event: Select.Changed) -> None:
        if event.value is Select.NULL:
            return
        tag = str(event.value)
        if tag not in self._mbc_available:
            return
        self._mbc_available.remove(tag)
        idx = self._next_mbc
        self._next_mbc += 1
        add_row = self.query_one("#mbc-add-row")
        self.query_one("#mechanical-bc-group").mount(MechanicalBCWidget(idx, tag, "mbc"), before=add_row)
        self.query_one("#mbc-add-select", Select).value = Select.NULL
        self._refresh_selects()

    def on_tag_released(self, event: TagReleased) -> None:
        if event.pool == "tbc":
            self._tbc_available.append(event.tag)
        elif event.pool == "mbc":
            self._mbc_available.append(event.tag)
        self._refresh_selects()

    def collect(self, problem: str) -> dict:
        tbc = {w._tag: w.collect() for w in self.query("#thermal-bc-group ThermalBCWidget").results(ThermalBCWidget)}
        mbc = {w._tag: w.collect() for w in self.query("#mechanical-bc-group MechanicalBCWidget").results(MechanicalBCWidget)}
        if problem == "thermal":
            return {"boundaries": tbc}
        if problem == "mechanical":
            return {"boundaries": mbc}
        return {"thermal_boundaries": tbc, "mechanical_boundaries": mbc}

    async def load(self, cfg: dict, problem: str) -> None:
        for w in list(self.query(ThermalBCWidget)):
            await w.remove()
        for w in list(self.query(MechanicalBCWidget)):
            await w.remove()

        tbc_data: dict = {}
        mbc_data: dict = {}
        if problem == "thermal":
            tbc_data = cfg.get("boundaries", {})
        elif problem == "mechanical":
            mbc_data = cfg.get("boundaries", {})
        else:
            tbc_data = cfg.get("thermal_boundaries", {})
            mbc_data = cfg.get("mechanical_boundaries", {})

        add_row_tbc = self.query_one("#tbc-add-row")
        for tag, data in tbc_data.items():
            if tag in self._tbc_available:
                self._tbc_available.remove(tag)
            idx = self._next_tbc
            self._next_tbc += 1
            w = ThermalBCWidget(idx, tag, "tbc")
            await self.query_one("#thermal-bc-group").mount(w, before=add_row_tbc)
            w.load(data)

        add_row_mbc = self.query_one("#mbc-add-row")
        for tag, data in mbc_data.items():
            if tag in self._mbc_available:
                self._mbc_available.remove(tag)
            idx = self._next_mbc
            self._next_mbc += 1
            w = MechanicalBCWidget(idx, tag, "mbc")
            await self.query_one("#mechanical-bc-group").mount(w, before=add_row_mbc)
            w.load(data)

        self._refresh_selects()

    def show_problem(self, problem: str) -> None:
        self._problem = problem
        self._apply_problem()

    def _apply_problem(self) -> None:
        tg = self.query_one("#thermal-bc-group")
        mg = self.query_one("#mechanical-bc-group")
        tl = self.query_one("#thermal-bc-label", Label)
        ml = self.query_one("#mechanical-bc-label", Label)
        p = self._problem
        if p == "thermal":
            tg.display = True
            mg.display = False
            tl.update("[bold][#5f87ff]boundaries[/][/]:")
        elif p == "mechanical":
            tg.display = False
            mg.display = True
            ml.update("[bold][#5f87ff]boundaries[/][/]:")
        else:
            tg.display = True
            mg.display = True
            tl.update("[bold][#5f87ff]thermal_boundaries[/][/]:")
            ml.update("[bold][#5f87ff]mechanical_boundaries[/][/]:")
