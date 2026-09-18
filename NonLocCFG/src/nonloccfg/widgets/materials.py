from __future__ import annotations

from textual import on
from textual.app import ComposeResult
from textual.containers import Horizontal, Vertical
from textual.widget import Widget
from textual.widgets import Button, Checkbox, Collapsible, Input, Label, RadioButton, RadioSet, Select

from nonloccfg.config_helpers import _format_value, _input, _parse_value, _radio, _select, _set_radio
from nonloccfg.messages import TagReleased
from nonloccfg.materials_db import MaterialsDB


class ModelFields(Widget):
    DEFAULT_CSS = """
    ModelFields { height: auto; }
    ModelFields .row { height: 1; align: left middle; margin-top: 1; }
    ModelFields .key { color: $text-muted; width: 20; }
    ModelFields .row Input, ModelFields .row Select { width: 1fr; }
    """

    def __init__(self, prefix: str) -> None:
        self._p = prefix
        super().__init__()

    def compose(self) -> ComposeResult:
        p = self._p
        with Horizontal(classes="row"):
            yield Label("local_weight:", classes="key")
            yield Input(placeholder="e.g. 0.78", id=f"{p}-lw", compact=True)
        with Horizontal(classes="row"):
            yield Label("nonlocal_radius:", classes="key")
            yield Input(placeholder="float or [rx, ry]", id=f"{p}-nlr", compact=True)
        with Horizontal(classes="row"):
            yield Label("search_radius:", classes="key")
            yield Input(placeholder="= nonlocal_radius", id=f"{p}-sr", compact=True)
        with Horizontal(classes="row"):
            yield Label("distance:", classes="key")
            with RadioSet(id=f"{p}-distance"):
                yield RadioButton("lp", value=True)
                yield RadioButton("ellipse_with_rotation")
        with Horizontal(classes="row"):
            yield Label("influence:", classes="key")
            yield Select(
                [("polynomial", "polynomial"), ("exponential", "exponential"), ("constant", "constant")],
                value="polynomial", id=f"{p}-influence", compact=True,
            )
        with Horizontal(classes="row"):
            yield Label("n:", classes="key")
            yield Input(placeholder="2  # default", id=f"{p}-n", compact=True)
        with Horizontal(classes="row"):
            yield Label("p:", classes="key")
            yield Input(placeholder="2  # default", id=f"{p}-p", compact=True)
        with Horizontal(classes="row"):
            yield Label("q:", classes="key")
            yield Input(placeholder="1.0  # default", id=f"{p}-q", compact=True)

    def collect(self) -> dict:
        p = self._p
        d: dict = {}
        for key, wid in (
            ("local_weight",    f"{p}-lw"),
            ("nonlocal_radius", f"{p}-nlr"),
            ("search_radius",   f"{p}-sr"),
            ("n", f"{p}-n"),
            ("p", f"{p}-p"),
            ("q", f"{p}-q"),
        ):
            v = _parse_value(_input(self, wid))
            if v is not None:
                d[key] = v
        distance = _radio(self, f"{p}-distance")
        if distance:
            d["distance"] = distance
        influence = _select(self, f"{p}-influence")
        if influence:
            d["influence"] = influence
        return d

    def load(self, data: dict) -> None:
        p = self._p
        for key, wid in (
            ("local_weight",    f"{p}-lw"),
            ("nonlocal_radius", f"{p}-nlr"),
            ("search_radius",   f"{p}-sr"),
            ("n", f"{p}-n"),
            ("p", f"{p}-p"),
            ("q", f"{p}-q"),
        ):
            if key in data:
                self.query_one(f"#{wid}", Input).value = _format_value(data[key])
        if "distance" in data:
            _set_radio(self, f"{p}-distance", data["distance"])
        if "influence" in data:
            self.query_one(f"#{p}-influence", Select).value = data["influence"]

    def set_locked(self, locked: bool) -> None:
        p = self._p
        for wid in (f"{p}-lw", f"{p}-nlr", f"{p}-sr", f"{p}-n", f"{p}-p", f"{p}-q"):
            try:
                self.query_one(f"#{wid}", Input).disabled = locked
            except Exception:
                pass
        # RadioSet and Select don't have disabled on the widget itself in all versions,
        # but Input fields cover the critical numeric ones.


class MaterialWidget(Widget):
    DEFAULT_CSS = """
    MaterialWidget { height: auto; padding-left: 2; margin-top: 1; }
    MaterialWidget .name-row { height: 1; align: left middle; }
    MaterialWidget .name-row Label { color: $text-muted; width: 2; }
    MaterialWidget .name-row Label.mat-tag { color: $text; width: auto; }
    MaterialWidget .name-row Button { width: 3; min-width: 3; }
    MaterialWidget .name-row Button.db-btn { width: 4; min-width: 4; margin-right: 1; }
    MaterialWidget .ref-row { height: 1; align: left middle; padding-left: 2; margin-top: 1; }
    MaterialWidget .ref-row Label { color: $accent; }
    MaterialWidget .ref-row Button { width: 9; min-width: 9; margin-left: 1; }
    MaterialWidget .sub-key { color: $text-muted; padding-left: 2; margin-top: 1; }
    MaterialWidget .row { height: 1; align: left middle; padding-left: 4; margin-top: 1; }
    MaterialWidget .key { color: $text-muted; width: 20; }
    MaterialWidget .row Input, MaterialWidget .row Select { width: 1fr; }
    MaterialWidget Vertical { height: auto; }
    MaterialWidget Collapsible { padding-left: 2; margin-top: 1; }
    MaterialWidget Checkbox { height: 1; margin-bottom: 1; padding: 0; }
    """

    def __init__(self, index: int, tag: str) -> None:
        self._index = index
        self._tag = tag
        self._db_ref: str | None = None
        super().__init__(id=f"mat-{index}")

    def compose(self) -> ComposeResult:
        i = self._index
        with Horizontal(classes="name-row"):
            yield Label("- ")
            yield Label(self._tag, classes="mat-tag")
            yield Button("DB", id=f"mat-db-{i}", compact=True, classes="db-btn")
            yield Button("×", variant="error", id=f"mat-remove-{i}", compact=True)
        with Horizontal(classes="ref-row", id=f"mat-ref-row-{i}"):
            yield Label("", id=f"mat-ref-label-{i}")
            yield Button("unlink", id=f"mat-unlink-{i}", variant="warning", compact=True)

        yield Label("physical:", classes="sub-key")

        with Vertical(id=f"mat-thermal-phys-{i}"):
            with Horizontal(classes="row"):
                yield Label("conductivity:", classes="key")
                yield Input(placeholder="float, [kx,ky] or 'x y: expr'", id=f"mat-conductivity-{i}", compact=True)
            with Horizontal(classes="row"):
                yield Label("capacity:", classes="key")
                yield Input(placeholder="1.0  # default", id=f"mat-capacity-{i}", compact=True)
            with Horizontal(classes="row"):
                yield Label("density:", classes="key")
                yield Input(placeholder="1.0  # default", id=f"mat-density-{i}", compact=True)
            with Horizontal(classes="row"):
                yield Label("relaxation_time:", classes="key")
                yield Input(placeholder="0.0  # default", id=f"mat-relax-{i}", compact=True)

        with Vertical(id=f"mat-mech-phys-{i}"):
            with Horizontal(classes="row"):
                yield Label("youngs_modulus:", classes="key")
                yield Input(placeholder="float or [Ex, Ey]", id=f"mat-youngs-{i}", compact=True)
            with Horizontal(classes="row"):
                yield Label("poissons_ratio:", classes="key")
                yield Input(placeholder="float or [null, v]", id=f"mat-poissons-{i}", compact=True)
            with Horizontal(classes="row"):
                yield Label("shear_modulus:", classes="key")
                yield Input(placeholder="float  # orthotropic only", id=f"mat-shear-{i}", compact=True)
            with Horizontal(classes="row"):
                yield Label("thermal_expansion:", classes="key")
                yield Input(placeholder="0.0  # default", id=f"mat-thexp-{i}", compact=True)

        with Collapsible(title="model:  # nonlocal, optional", collapsed=True, id=f"mat-model-{i}"):
            yield ModelFields(prefix=f"mat-{i}-m")
        with Collapsible(title="thermal_model:  # nonlocal, optional", collapsed=True, id=f"mat-thermal-model-{i}"):
            yield ModelFields(prefix=f"mat-{i}-tm")
        with Collapsible(title="mechanical_model:  # nonlocal, optional", collapsed=True, id=f"mat-mech-model-{i}"):
            yield Checkbox("same as thermal_model", id=f"mat-mech-same-{i}", value=False, compact=True)
            yield ModelFields(prefix=f"mat-{i}-mm")

    def on_mount(self) -> None:
        self.query_one(f"#mat-ref-row-{self._index}").display = False
        try:
            sel = self.app.query_one("#task-problem", Select)
            problem = "thermal" if sel.value is Select.NULL else str(sel.value)
        except Exception:
            problem = "thermal"
        self.show_problem(problem)

    def show_problem(self, problem: str) -> None:
        i = self._index
        is_thermo = problem == "thermomechanical"
        self.query_one(f"#mat-thermal-phys-{i}").display = problem in ("thermal", "thermomechanical")
        self.query_one(f"#mat-mech-phys-{i}").display = problem in ("mechanical", "thermomechanical")
        self.query_one(f"#mat-model-{i}").display = not is_thermo
        self.query_one(f"#mat-thermal-model-{i}").display = is_thermo
        self.query_one(f"#mat-mech-model-{i}").display = is_thermo

    def on_checkbox_changed(self, event: Checkbox.Changed) -> None:
        if event.checkbox.id == f"mat-mech-same-{self._index}":
            self._apply_mech_same(event.value)

    def _apply_mech_same(self, same: bool) -> None:
        i = self._index
        try:
            self.query_one(f"#mat-mech-model-{i}", Collapsible).query_one(ModelFields).display = not same
        except Exception:
            pass

    def on_button_pressed(self, event: Button.Pressed) -> None:
        i = self._index
        if event.button.id == f"mat-remove-{i}":
            event.stop()
            self.post_message(TagReleased(self._tag, "mat"))
            self.remove()
        elif event.button.id == f"mat-db-{i}":
            event.stop()
            from nonloccfg.widgets.material_picker import MaterialPickerModal
            self.app.push_screen(MaterialPickerModal(), self._on_db_picked)
        elif event.button.id == f"mat-unlink-{i}":
            event.stop()
            self._unlink()

    def _on_db_picked(self, name: str | None) -> None:
        if not name:
            return
        from nonloccfg.db import db
        data = db.get(name)
        if not data:
            return
        self._db_ref = name
        i = self._index
        # Fill and lock physical fields
        for field in MaterialsDB.FIELDS:
            wid_id = f"mat-{self._FIELD_TO_WID.get(field, field)}-{i}"
            try:
                inp = self.query_one(f"#{wid_id}", Input)
                inp.value = _format_value(data[field]) if field in data else ""
                inp.disabled = True
            except Exception:
                pass
        # Fill and lock model fields based on problem type
        try:
            sel = self.app.query_one("#task-problem", Select)
            problem = "thermal" if sel.value is Select.NULL else str(sel.value)
        except Exception:
            problem = "thermal"

        if problem == "thermomechanical":
            tm_data = data.get("thermal_model")
            mm_data = data.get("mechanical_model") or tm_data
            for coll_id, mdata in (
                (f"mat-thermal-model-{i}", tm_data),
                (f"mat-mech-model-{i}",    mm_data),
            ):
                if mdata:
                    try:
                        coll = self.query_one(f"#{coll_id}", Collapsible)
                        coll.collapsed = False
                        mf = coll.query_one(ModelFields)
                        mf.load(mdata)
                        mf.set_locked(True)
                    except Exception:
                        pass
        else:
            db_key = "mechanical_model" if problem == "mechanical" else "thermal_model"
            model_data = data.get(db_key)
            if model_data:
                try:
                    coll = self.query_one(f"#mat-model-{i}", Collapsible)
                    coll.collapsed = False
                    mf = coll.query_one(ModelFields)
                    mf.load(model_data)
                    mf.set_locked(True)
                except Exception:
                    pass
        # Show ref row
        self.query_one(f"#mat-ref-label-{i}", Label).update(f"ref: [bold]{name}[/]")
        self.query_one(f"#mat-ref-row-{i}").display = True

    def _unlink(self) -> None:
        self._db_ref = None
        i = self._index
        for field in MaterialsDB.FIELDS:
            wid_id = f"mat-{self._FIELD_TO_WID.get(field, field)}-{i}"
            try:
                self.query_one(f"#{wid_id}", Input).disabled = False
            except Exception:
                pass
        for coll_id in (f"mat-model-{i}", f"mat-thermal-model-{i}", f"mat-mech-model-{i}"):
            try:
                self.query_one(f"#{coll_id}", Collapsible).query_one(ModelFields).set_locked(False)
            except Exception:
                pass
        self.query_one(f"#mat-ref-row-{i}").display = False

    # Mapping from DB field names to widget id fragments (where they differ)
    _FIELD_TO_WID: dict[str, str] = {
        "conductivity":    "conductivity",
        "capacity":        "capacity",
        "density":         "density",
        "relaxation_time": "relax",
        "youngs_modulus":  "youngs",
        "poissons_ratio":  "poissons",
        "shear_modulus":   "shear",
        "thermal_expansion": "thexp",
    }

    def load(self, data: dict, problem: str) -> None:
        i = self._index
        # Restore DB ref state if present
        if "$db_ref" in data:
            self._db_ref = data["$db_ref"]
            self.query_one(f"#mat-ref-label-{i}", Label).update(
                f"ref: [bold]{self._db_ref}[/]"
            )
            self.query_one(f"#mat-ref-row-{i}").display = True
        phys = data.get("physical", {})
        for key, wid in (
            ("conductivity",      f"mat-conductivity-{i}"),
            ("capacity",          f"mat-capacity-{i}"),
            ("density",           f"mat-density-{i}"),
            ("relaxation_time",   f"mat-relax-{i}"),
            ("youngs_modulus",    f"mat-youngs-{i}"),
            ("poissons_ratio",    f"mat-poissons-{i}"),
            ("shear_modulus",     f"mat-shear-{i}"),
            ("thermal_expansion", f"mat-thexp-{i}"),
        ):
            if key in phys:
                inp = self.query_one(f"#{wid}", Input)
                inp.value = _format_value(phys[key])
                if self._db_ref:
                    inp.disabled = True

        if problem == "thermomechanical":
            if "thermal_model" in data:
                coll = self.query_one(f"#mat-thermal-model-{i}", Collapsible)
                coll.collapsed = False
                mf = coll.query_one(ModelFields)
                mf.load(data["thermal_model"])
                if self._db_ref:
                    mf.set_locked(True)
            mm_coll = self.query_one(f"#mat-mech-model-{i}", Collapsible)
            if "mechanical_model" in data:
                same = data.get("mechanical_model") == data.get("thermal_model")
                cb = self.query_one(f"#mat-mech-same-{i}", Checkbox)
                cb.value = same
                self._apply_mech_same(same)
                if not same:
                    mm_coll.collapsed = False
                    mf = mm_coll.query_one(ModelFields)
                    mf.load(data["mechanical_model"])
                    if self._db_ref:
                        mf.set_locked(True)
        else:
            if "model" in data:
                coll = self.query_one(f"#mat-model-{i}", Collapsible)
                coll.collapsed = False
                mf = coll.query_one(ModelFields)
                mf.load(data["model"])
                if self._db_ref:
                    mf.set_locked(True)

    def collect(self, problem: str) -> dict:
        i = self._index
        d: dict = {}
        if self._db_ref:
            d["$db_ref"] = self._db_ref

        phys: dict = {}
        if problem in ("thermal", "thermomechanical"):
            for key, wid in (
                ("conductivity",    f"mat-conductivity-{i}"),
                ("capacity",        f"mat-capacity-{i}"),
                ("density",         f"mat-density-{i}"),
                ("relaxation_time", f"mat-relax-{i}"),
            ):
                v = _parse_value(_input(self, wid))
                if v is not None:
                    phys[key] = v
        if problem in ("mechanical", "thermomechanical"):
            for key, wid in (
                ("youngs_modulus",    f"mat-youngs-{i}"),
                ("poissons_ratio",    f"mat-poissons-{i}"),
                ("shear_modulus",     f"mat-shear-{i}"),
                ("thermal_expansion", f"mat-thexp-{i}"),
            ):
                v = _parse_value(_input(self, wid))
                if v is not None:
                    phys[key] = v
        if phys:
            d["physical"] = phys

        if problem == "thermomechanical":
            tm = self.query_one(f"#mat-thermal-model-{i}", Collapsible)
            mm = self.query_one(f"#mat-mech-model-{i}", Collapsible)
            tm_data = tm.query_one(ModelFields).collect()
            if "local_weight" in tm_data and "nonlocal_radius" in tm_data:
                d["thermal_model"] = tm_data
            same = self.query_one(f"#mat-mech-same-{i}", Checkbox).value
            if same:
                if "thermal_model" in d:
                    d["mechanical_model"] = d["thermal_model"]
            else:
                mm_data = mm.query_one(ModelFields).collect()
                if "local_weight" in mm_data and "nonlocal_radius" in mm_data:
                    d["mechanical_model"] = mm_data
        else:
            m = self.query_one(f"#mat-model-{i}", Collapsible)
            m_data = m.query_one(ModelFields).collect()
            if "local_weight" in m_data and "nonlocal_radius" in m_data:
                d["model"] = m_data

        return d


class MaterialsSection(Widget):
    DEFAULT_CSS = """
    MaterialsSection { height: auto; margin-bottom: 1; }
    MaterialsSection .add-row { height: 1; align: left middle; padding-left: 2; margin-top: 1; }
    MaterialsSection .add-row Label { color: $text-muted; width: 6; }
    MaterialsSection .add-row Select { width: 1fr; }
    """

    def __init__(self) -> None:
        super().__init__()
        self._available: list[str] = []
        self._next_index = 0

    def compose(self) -> ComposeResult:
        yield Label("[bold][#5f87ff]materials[/][/]:")
        with Horizontal(classes="add-row", id="mat-add-row"):
            yield Label("add:")
            yield Select([], allow_blank=True, prompt="— select material —",
                         id="mat-add-select", compact=True)

    def on_mount(self) -> None:
        self.query_one("#mat-add-select", Select).disabled = True

    def set_material_tags(self, tags: list[str]) -> None:
        self._available = list(tags)
        self._refresh_select()

    def _refresh_select(self) -> None:
        sel = self.query_one("#mat-add-select", Select)
        sel.set_options([(t, t) for t in self._available])
        sel.disabled = len(self._available) == 0

    @on(Select.Changed, "#mat-add-select")
    def _add_mat(self, event: Select.Changed) -> None:
        if event.value is Select.NULL:
            return
        tag = str(event.value)
        if tag not in self._available:
            return
        self._available.remove(tag)
        idx = self._next_index
        self._next_index += 1
        add_row = self.query_one("#mat-add-row")
        self.mount(MaterialWidget(idx, tag), before=add_row)
        self.query_one("#mat-add-select", Select).value = Select.NULL
        self._refresh_select()

    def on_tag_released(self, event: TagReleased) -> None:
        if event.pool == "mat":
            self._available.append(event.tag)
            self._refresh_select()

    def collect(self, problem: str) -> dict:
        return {w._tag: w.collect(problem) for w in self.query(MaterialWidget)}

    async def load(self, data: dict, problem: str) -> None:
        for w in list(self.query(MaterialWidget)):
            await w.remove()
        add_row = self.query_one("#mat-add-row")
        for tag, mat_data in data.items():
            if tag in self._available:
                self._available.remove(tag)
            idx = self._next_index
            self._next_index += 1
            w = MaterialWidget(idx, tag)
            await self.mount(w, before=add_row)
            w.load(mat_data, problem)
        self._refresh_select()

    def show_problem(self, problem: str) -> None:
        for w in self.query(MaterialWidget):
            w.show_problem(problem)
