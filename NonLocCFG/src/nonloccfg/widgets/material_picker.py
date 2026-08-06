from __future__ import annotations

from textual import on
from textual.app import ComposeResult
from textual.containers import Horizontal, Vertical, VerticalScroll
from textual.screen import ModalScreen
from textual.widgets import Button, Input, Label, ListItem, ListView

from nonloccfg.db import db
from nonloccfg.materials_db import MaterialsDB


class MaterialPickerModal(ModalScreen[str | None]):
    DEFAULT_CSS = """
    MaterialPickerModal { align: center middle; }
    MaterialPickerModal #dialog {
        width: 70; height: 28;
        border: solid $accent;
        background: $surface;
        padding: 1;
    }
    MaterialPickerModal #search-row {
        height: 1; margin-bottom: 1;
    }
    MaterialPickerModal #search-row Label {
        width: 10; color: $text-muted;
    }
    MaterialPickerModal #search-row Input { width: 1fr; }
    MaterialPickerModal #body { height: 1fr; }
    MaterialPickerModal #list-col { width: 26; }
    MaterialPickerModal #preview-col { width: 1fr; padding-left: 1; }
    MaterialPickerModal ListView { height: 1fr; border: solid $panel; }
    MaterialPickerModal #preview-scroll { height: 1fr; border: solid $panel; padding: 0 1; }
    MaterialPickerModal .preview-key { color: $text-muted; width: 22; }
    MaterialPickerModal .preview-row { height: 1; }
    MaterialPickerModal #btn-row {
        height: 1; margin-top: 1; align: right middle;
    }
    MaterialPickerModal #btn-row Button { margin-left: 1; }
    """
    BINDINGS = [("escape", "dismiss(None)", "Cancel")]

    def compose(self) -> ComposeResult:
        with Vertical(id="dialog"):
            with Horizontal(id="search-row"):
                yield Label("search:", )
                yield Input(placeholder="type to filter…", id="search-input", compact=True)
            with Horizontal(id="body"):
                with Vertical(id="list-col"):
                    yield ListView(id="mat-list")
                with Vertical(id="preview-col"):
                    yield Label("[dim]select a material[/]", id="preview-title")
                    with VerticalScroll(id="preview-scroll"):
                        pass
            with Horizontal(id="btn-row"):
                yield Button("Cancel", id="btn-cancel", variant="default", compact=True)
                yield Button("Select", id="btn-select", variant="primary", disabled=True, compact=True)

    def on_mount(self) -> None:
        self._rebuild_list(db.all_names())

    def _rebuild_list(self, names: list[str]) -> None:
        lv = self.query_one("#mat-list", ListView)
        lv.clear()
        for name in names:
            lv.append(ListItem(Label(name), name=name))

    @on(Input.Changed, "#search-input")
    def _on_search(self, event: Input.Changed) -> None:
        q = event.value.strip()
        results = db.search(q) if q else db.all_names()
        self._rebuild_list(results)
        self._clear_preview()
        self.query_one("#btn-select", Button).disabled = True

    @on(ListView.Highlighted)
    def _on_highlighted(self, event: ListView.Highlighted) -> None:
        if event.item is None:
            self._clear_preview()
            return
        name = event.item.name
        self._show_preview(name)
        self.query_one("#btn-select", Button).disabled = False

    def _clear_preview(self) -> None:
        self.query_one("#preview-title", Label).update("[dim]select a material[/]")
        scroll = self.query_one("#preview-scroll", VerticalScroll)
        for child in list(scroll.children):
            child.remove()

    def _show_preview(self, name: str) -> None:
        data = db.get(name)
        self.query_one("#preview-title", Label).update(f"[bold]{name}[/]")
        scroll = self.query_one("#preview-scroll", VerticalScroll)
        for child in list(scroll.children):
            child.remove()
        if not data:
            return
        for field in MaterialsDB.FIELDS:
            if field in data:
                row = Horizontal(classes="preview-row")
                scroll.mount(row)
                row.mount(Label(f"{field}:", classes="preview-key"))
                row.mount(Label(str(data[field])))
        for mkey in ("thermal_model", "mechanical_model"):
            model = data.get(mkey)
            if model:
                scroll.mount(Label(f"[dim]{mkey}:[/]"))
                for field in MaterialsDB.MODEL_FIELDS:
                    if field in model:
                        row = Horizontal(classes="preview-row")
                        scroll.mount(row)
                        row.mount(Label(f"  {field}:", classes="preview-key"))
                        row.mount(Label(str(model[field])))

    @on(Button.Pressed, "#btn-select")
    def _select(self) -> None:
        lv = self.query_one("#mat-list", ListView)
        if lv.highlighted_child is not None:
            self.dismiss(lv.highlighted_child.name)
        else:
            self.dismiss(None)

    @on(Button.Pressed, "#btn-cancel")
    def _cancel(self) -> None:
        self.dismiss(None)
