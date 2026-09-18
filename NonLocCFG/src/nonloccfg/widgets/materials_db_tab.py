from __future__ import annotations

from textual import on
from textual.app import ComposeResult
from textual.containers import Horizontal, Vertical, VerticalScroll
from textual.widget import Widget
from textual.widgets import Button, Checkbox, Input, Label, ListItem, ListView, RadioButton, RadioSet, Select

from nonloccfg.config_helpers import _format_value, _parse_value, _radio, _select, _set_radio
from nonloccfg.db import db
from nonloccfg.materials_db import MaterialsDB


class MaterialsDBTab(Widget):
    DEFAULT_CSS = """
    MaterialsDBTab { height: 1fr; }
    MaterialsDBTab #layout { height: 1fr; }

    MaterialsDBTab #left-col {
        width: 28; border-right: solid $panel; padding: 0 1; margin-bottom: 1;
    }
    MaterialsDBTab #left-col Label.section-title {
        color: $text-muted; margin-top: 1; margin-bottom: 1;
    }
    MaterialsDBTab ListView { height: 1fr; }
    MaterialsDBTab #left-buttons {
        height: 1; margin-top: 1;
    }
    MaterialsDBTab #btn-new { width: 1fr; min-width: 1; }
    MaterialsDBTab #btn-delete { width: 1fr; min-width: 1; margin-left: 1; }

    MaterialsDBTab #right-col { width: 1fr; padding: 0 2; }
    MaterialsDBTab #editor-title {
        margin-top: 1; margin-bottom: 1;
    }
    MaterialsDBTab #name-row {
        height: 1; align: left middle; margin-bottom: 1;
    }
    MaterialsDBTab #name-row Label { color: $text-muted; width: 22; }
    MaterialsDBTab #name-row Input { width: 1fr; }
    MaterialsDBTab .field-row {
        height: 1; align: left middle; margin-top: 1;
    }
    MaterialsDBTab .field-key { color: $text-muted; width: 22; }
    MaterialsDBTab .field-row Input { width: 1fr; }
    MaterialsDBTab .section-label {
        color: $text-muted; margin-top: 1; padding-left: 0;
    }
    MaterialsDBTab Checkbox { height: 1; margin-top: 1; padding: 0; }
    MaterialsDBTab .field-row RadioSet { width: 1fr; }
    MaterialsDBTab .field-row Select { width: 1fr; }
    MaterialsDBTab #editor-buttons {
        height: 1; margin-top: 1; align: right middle;
    }
    MaterialsDBTab #editor-buttons Button { margin-left: 1; }
    MaterialsDBTab #no-selection {
        color: $text-muted; margin-top: 2; padding-left: 2;
    }
    """

    def compose(self) -> ComposeResult:
        with Horizontal(id="layout"):
            with Vertical(id="left-col"):
                yield Label("Materials DB", classes="section-title")
                yield ListView(id="db-list")
                with Horizontal(id="left-buttons"):
                    yield Button("New", id="btn-new", compact=True)
                    yield Button("Delete", id="btn-delete", variant="error", compact=True, disabled=True)
            with Vertical(id="right-col"):
                yield Label("[dim]Select a material or press New[/]", id="no-selection")
                with VerticalScroll(id="editor"):
                    yield Label("", id="editor-title")
                    with Horizontal(id="name-row"):
                        yield Label("name:", classes="field-key")
                        yield Input(placeholder="Material name", id="field-name", compact=True)
                    yield Label("[bold][#5f87ff]physical[/][/]", classes="section-label")
                    for field in MaterialsDB.FIELDS:
                        with Horizontal(classes="field-row", id=f"frow-{field}"):
                            yield Label(f"{field}:", classes="field-key")
                            yield Input(placeholder="float  # optional", id=f"field-{field}", compact=True)
                    for prefix, title in (
                        ("thermal_model",    "[bold][#5f87ff]thermal_model[/][/]"),
                        ("mechanical_model", "[bold][#5f87ff]mechanical_model[/][/]"),
                    ):
                        yield Label(title, classes="section-label")
                        if prefix == "mechanical_model":
                            Checkbox.BUTTON_INNER = "v"
                            yield Checkbox("same as thermal_model", id="mech-same-as-thermal", value=True, compact=True)
                        row_cls = f"field-row mrow-{prefix}"
                        with Horizontal(classes=row_cls):
                            yield Label("local_weight:", classes="field-key")
                            yield Input(placeholder="e.g. 0.78", id=f"field-{prefix}-lw", compact=True)
                        with Horizontal(classes=row_cls):
                            yield Label("nonlocal_radius:", classes="field-key")
                            yield Input(placeholder="float or [rx, ry]", id=f"field-{prefix}-nlr", compact=True)
                        with Horizontal(classes=row_cls):
                            yield Label("search_radius:", classes="field-key")
                            yield Input(placeholder="= nonlocal_radius", id=f"field-{prefix}-sr", compact=True)
                        with Horizontal(classes=row_cls):
                            yield Label("distance:", classes="field-key")
                            with RadioSet(id=f"field-{prefix}-distance"):
                                yield RadioButton("lp", value=True)
                                yield RadioButton("ellipse_with_rotation")
                        with Horizontal(classes=row_cls):
                            yield Label("influence:", classes="field-key")
                            yield Select(
                                [("polynomial", "polynomial"), ("exponential", "exponential"), ("constant", "constant")],
                                value="polynomial", id=f"field-{prefix}-influence", compact=True,
                            )
                        with Horizontal(classes=row_cls):
                            yield Label("n:", classes="field-key")
                            yield Input(placeholder="2  # default", id=f"field-{prefix}-n", compact=True)
                        with Horizontal(classes=row_cls):
                            yield Label("p:", classes="field-key")
                            yield Input(placeholder="2  # default", id=f"field-{prefix}-p", compact=True)
                        with Horizontal(classes=row_cls):
                            yield Label("q:", classes="field-key")
                            yield Input(placeholder="1.0  # default", id=f"field-{prefix}-q", compact=True)
                    with Horizontal(id="editor-buttons"):
                        yield Button("Save", id="btn-save", variant="primary", compact=True)

    def on_mount(self) -> None:
        self.query_one("#editor").display = False
        self._refresh_list()
        self._apply_same_as_thermal(True)

    @on(Checkbox.Changed, "#mech-same-as-thermal")
    def _on_same_toggled(self, event: Checkbox.Changed) -> None:
        self._apply_same_as_thermal(event.value)

    def _apply_same_as_thermal(self, same: bool) -> None:
        for row in self.query(".mrow-mechanical_model"):
            row.display = not same

    def _refresh_list(self, select_name: str | None = None) -> None:
        lv = self.query_one("#db-list", ListView)
        lv.clear()
        names = db.all_names()
        for name in names:
            lv.append(ListItem(Label(name), name=name))
        if select_name and select_name in names:
            lv.index = names.index(select_name)

    def _show_editor(self, name: str | None) -> None:
        self.query_one("#no-selection").display = False
        editor = self.query_one("#editor")
        was_hidden = not editor.display
        editor.display = True
        if was_hidden:
            # Force layout recalculation after making visible
            self.call_after_refresh(editor.refresh, layout=True)

        if name is None:
            self.query_one("#editor-title", Label).update("[bold]New material[/]")
            self.query_one("#field-name", Input).value = ""
            for field in MaterialsDB.FIELDS:
                self.query_one(f"#field-{field}", Input).value = ""
            for prefix in ("thermal_model", "mechanical_model"):
                self._clear_model_section(prefix)
            self.query_one("#mech-same-as-thermal", Checkbox).value = True
            self._apply_same_as_thermal(True)
        else:
            self.query_one("#editor-title", Label).update(f"[bold]{name}[/]")
            self.query_one("#field-name", Input).value = name
            data = db.get(name) or {}
            for field in MaterialsDB.FIELDS:
                v = data.get(field)
                self.query_one(f"#field-{field}", Input).value = _format_value(v) if v is not None else ""
            for prefix in ("thermal_model", "mechanical_model"):
                self._load_model_section(prefix, data.get(prefix) or {})
            same = "mechanical_model" not in data or data.get("mechanical_model") == data.get("thermal_model")
            self.query_one("#mech-same-as-thermal", Checkbox).value = same
            self._apply_same_as_thermal(same)

    @on(ListView.Highlighted)
    def _on_highlighted(self, event: ListView.Highlighted) -> None:
        has = event.item is not None
        self.query_one("#btn-delete", Button).disabled = not has
        if has:
            self._show_editor(event.item.name)

    @on(Button.Pressed, "#btn-new")
    def _new(self) -> None:
        self.query_one("#db-list", ListView).index = None
        self._show_editor(None)
        self.query_one("#btn-delete", Button).disabled = True
        self.query_one("#field-name", Input).focus()

    @on(Button.Pressed, "#btn-delete")
    def _delete(self) -> None:
        lv = self.query_one("#db-list", ListView)
        if lv.highlighted_child is not None:
            db.delete(lv.highlighted_child.name)
            self._refresh_list()
            self.query_one("#no-selection").display = True
            self.query_one("#editor").display = False
            self.query_one("#btn-delete", Button).disabled = True

    def _clear_model_section(self, prefix: str) -> None:
        for wid in (f"field-{prefix}-lw", f"field-{prefix}-nlr", f"field-{prefix}-sr",
                    f"field-{prefix}-n", f"field-{prefix}-p", f"field-{prefix}-q"):
            self.query_one(f"#{wid}", Input).value = ""

    def _load_model_section(self, prefix: str, model: dict) -> None:
        for key, wid in (
            ("local_weight",    f"field-{prefix}-lw"),
            ("nonlocal_radius", f"field-{prefix}-nlr"),
            ("search_radius",   f"field-{prefix}-sr"),
            ("n", f"field-{prefix}-n"),
            ("p", f"field-{prefix}-p"),
            ("q", f"field-{prefix}-q"),
        ):
            v = model.get(key)
            self.query_one(f"#{wid}", Input).value = _format_value(v) if v is not None else ""
        if "distance" in model:
            _set_radio(self, f"field-{prefix}-distance", model["distance"])
        if "influence" in model:
            self.query_one(f"#field-{prefix}-influence", Select).value = model["influence"]

    def _collect_model_section(self, prefix: str) -> dict:
        model: dict = {}
        for key, wid in (
            ("local_weight",    f"field-{prefix}-lw"),
            ("nonlocal_radius", f"field-{prefix}-nlr"),
            ("search_radius",   f"field-{prefix}-sr"),
            ("n", f"field-{prefix}-n"),
            ("p", f"field-{prefix}-p"),
            ("q", f"field-{prefix}-q"),
        ):
            v = _parse_value(self.query_one(f"#{wid}", Input).value)
            if v is not None:
                model[key] = v
        distance = _radio(self, f"field-{prefix}-distance")
        if distance and distance != "lp":
            model["distance"] = distance
        influence = _select(self, f"field-{prefix}-influence")
        if influence and influence != "polynomial":
            model["influence"] = influence
        return model

    @on(Button.Pressed, "#btn-save")
    def _save(self) -> None:
        name = self.query_one("#field-name", Input).value.strip()
        if not name:
            self.app.notify("Name is required", severity="warning")
            return
        data: dict = {}
        for field in MaterialsDB.FIELDS:
            v = _parse_value(self.query_one(f"#field-{field}", Input).value)
            if v is not None:
                data[field] = v
        m = self._collect_model_section("thermal_model")
        if "local_weight" in m and "nonlocal_radius" in m:
            data["thermal_model"] = m
        same = self.query_one("#mech-same-as-thermal", Checkbox).value
        if not same:
            m = self._collect_model_section("mechanical_model")
            if "local_weight" in m and "nonlocal_radius" in m:
                data["mechanical_model"] = m
        db.upsert(name, data)
        self._refresh_list(select_name=name)
        self._show_editor(name)
        self.app.notify(f"Saved '{name}'", severity="information")
