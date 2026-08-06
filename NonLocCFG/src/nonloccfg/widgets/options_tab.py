from __future__ import annotations

from textual import on
from textual.app import ComposeResult
from textual.containers import Horizontal, VerticalScroll
from textual.widget import Widget
from textual.widgets import Button, Input, Label

from nonloccfg.config_helpers import _format_value, _input, _parse_value
from nonloccfg.widgets.file_picker import FilePickerModal


class SaveSection(Widget):
    DEFAULT_CSS = """
    SaveSection { height: auto; margin-bottom: 1; }
    SaveSection .row { height: 1; align: left middle; padding-left: 2; margin-top: 1; }
    SaveSection .key { color: $text-muted; width: 14; }
    SaveSection .row Input { width: 1fr; }
    """

    def compose(self) -> ComposeResult:
        yield Label("[bold][#5f87ff]save[/][/]:")
        for key, wid, ph in (
            ("folder:",    "save-folder",    "./results"),
            ("config:",    "save-config",    "config  # filename without ext"),
            ("csv:",       "save-csv",       "solution  # filename without ext"),
            ("vtk:",       "save-vtk",       "solution  # filename without ext"),
            ("precision:", "save-precision", "7  # default"),
        ):
            with Horizontal(classes="row"):
                yield Label(key, classes="key")
                yield Input(placeholder=ph, id=wid, compact=True)

    def collect(self) -> dict:
        d: dict = {}
        for key, wid in (
            ("folder",    "save-folder"),
            ("config",    "save-config"),
            ("csv",       "save-csv"),
            ("vtk",       "save-vtk"),
            ("precision", "save-precision"),
        ):
            v = _parse_value(_input(self, wid))
            if v is not None:
                d[key] = v
        return d

    def load(self, data: dict) -> None:
        for key, wid in (
            ("folder",    "save-folder"),
            ("config",    "save-config"),
            ("csv",       "save-csv"),
            ("vtk",       "save-vtk"),
            ("precision", "save-precision"),
        ):
            if key in data:
                self.query_one(f"#{wid}", Input).value = _format_value(data[key])


class AuxiliarySection(Widget):
    DEFAULT_CSS = """
    AuxiliarySection { height: auto; margin-bottom: 1; }
    AuxiliarySection .tbd { color: $text-muted; padding-left: 2; margin-top: 1; }
    """

    def compose(self) -> ComposeResult:
        yield Label("[bold][#5f87ff]auxiliary[/][/]:")
        yield Label("TBD", classes="tbd")

    def collect(self) -> dict:
        return {}

    def load(self, data: dict) -> None:
        pass


class OutputSection(Widget):
    DEFAULT_CSS = """
    OutputSection { height: auto; margin-bottom: 1; }
    OutputSection .row { height: 1; align: left middle; padding-left: 2; margin-top: 1; }
    OutputSection .key { color: $text-muted; width: 14; }
    OutputSection .row Input { width: 1fr; }
    OutputSection .row Button { width: 10; min-width: 10; }
    """

    def __init__(self) -> None:
        super().__init__()
        self._current_path: str = ""

    def compose(self) -> ComposeResult:
        yield Label("[bold][#5f87ff]config output[/][/]:")
        with Horizontal(classes="row"):
            yield Label("path:", classes="key")
            yield Input(placeholder="config.json", id="output-path", compact=True)
            yield Button("Browse", id="output-browse", compact=True)

    @on(Button.Pressed, "#output-browse")
    def _browse(self) -> None:
        self.app.push_screen(
            FilePickerModal(self._current_path, suffixes=frozenset({".json"})),
            self._on_picked,
        )

    def _on_picked(self, path: str | None) -> None:
        if path:
            self._current_path = path
            self.query_one("#output-path", Input).value = path

    def get_path(self) -> str:
        return _input(self, "output-path").strip() or "config.json"


class OptionsTab(Widget):
    DEFAULT_CSS = """
    OptionsTab { height: 1fr; }
    OptionsTab VerticalScroll { padding: 0 2; }
    """

    def compose(self) -> ComposeResult:
        with VerticalScroll():
            yield SaveSection()
            yield AuxiliarySection()
            yield OutputSection()

    def build_options(self) -> dict:
        return {
            "save":      self.query_one(SaveSection).collect(),
            "auxiliary": self.query_one(AuxiliarySection).collect(),
        }

    def load_options(self, cfg: dict) -> None:
        if "save" in cfg:
            self.query_one(SaveSection).load(cfg["save"])
        if "auxiliary" in cfg:
            self.query_one(AuxiliarySection).load(cfg["auxiliary"])

    def get_output_path(self) -> str:
        return self.query_one(OutputSection).get_path()
