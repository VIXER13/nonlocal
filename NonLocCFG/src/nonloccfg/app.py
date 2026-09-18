from __future__ import annotations

import json
from pathlib import Path

from textual.app import App, ComposeResult
from textual.widgets import Footer, Header, TabbedContent, TabPane

from nonloccfg.widgets.file_picker import FilePickerModal
from nonloccfg.widgets.materials_db_tab import MaterialsDBTab
from nonloccfg.widgets.options_tab import OptionsTab
from nonloccfg.widgets.task_tab import TaskTab


class NonlocCfgApp(App):
    CSS = """
    RadioSet {
        layout: horizontal;
        height: 1;
        border: none;
        padding: 0;
        margin: 0;
        background: transparent;
        width: auto;
    }
    RadioButton {
        height: 1;
        width: auto;
        padding: 0 1;
        background: transparent;
    }
    RadioButton > .toggle--label { width: auto; }
    """
    BINDINGS = [
        ("q", "quit", "Quit"),
        ("ctrl+o", "open_config", "Open config"),
        ("ctrl+s", "dump_config", "Save config"),
    ]

    def compose(self) -> ComposeResult:
        yield Header()
        with TabbedContent():
            with TabPane("Task", id="tab-task"):
                yield TaskTab()
            with TabPane("Materials DB", id="tab-materials"):
                yield MaterialsDBTab()
            with TabPane("Options", id="tab-options"):
                yield OptionsTab()
        yield Footer()

    def action_open_config(self) -> None:
        start = self.query_one(OptionsTab).get_output_path()
        self.push_screen(
            FilePickerModal(start_path=start, suffixes=frozenset({".json"})),
            self._on_config_picked,
        )

    def _on_config_picked(self, path: str | None) -> None:
        if not path:
            return
        config_path = Path(path).resolve()
        try:
            cfg = json.loads(config_path.read_text(encoding="utf-8"))
        except Exception as exc:
            self.notify(f"Read error: {exc}", severity="error")
            return
        self.query_one(OptionsTab).load_options(cfg)
        self.run_worker(
            self.query_one(TaskTab).load_config(cfg, config_path.parent),
            exclusive=True,
        )

    def action_dump_config(self) -> None:
        try:
            opts_tab = self.query_one(OptionsTab)
            cfg = self.query_one(TaskTab).build_config()
            opts = opts_tab.build_options()
            if opts["save"]:
                cfg["save"] = opts["save"]
            if opts["auxiliary"]:
                cfg["auxiliary"] = opts["auxiliary"]
            out_path = Path(opts_tab.get_output_path())
        except Exception as exc:
            self.notify(f"Build error: {exc}", severity="error")
            return
        try:
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(json.dumps(cfg, indent=4, ensure_ascii=False), encoding="utf-8")
            self.notify(f"Saved → {out_path.absolute()}", severity="information")
        except OSError as exc:
            self.notify(f"Write error: {exc}", severity="error")


def main() -> None:
    NonlocCfgApp().run()
