from __future__ import annotations

from pathlib import Path
from typing import Iterable

from textual import on
from textual.app import ComposeResult
from textual.containers import Horizontal, Vertical
from textual.screen import ModalScreen
from textual.widgets import Button, DirectoryTree, Input, Label


class FilteredTree(DirectoryTree):
    """DirectoryTree that shows only directories and files with given suffixes."""

    def __init__(self, path: str, suffixes: frozenset[str], **kwargs) -> None:
        self._suffixes = suffixes
        super().__init__(path, **kwargs)

    def filter_paths(self, paths: Iterable[Path]) -> Iterable[Path]:
        return [p for p in paths if p.is_dir() or p.suffix.lower() in self._suffixes]


class FilePickerModal(ModalScreen[str | None]):
    DEFAULT_CSS = """
    FilePickerModal {
        align: center middle;
    }
    FilePickerModal #dialog {
        width: 80;
        height: 30;
        border: solid $accent;
        background: $surface;
        padding: 1;
    }
    FilePickerModal #nav-row {
        height: 1;
        align: left middle;
        margin-bottom: 1;
    }
    FilePickerModal #nav-row Button {
        width: 4;
        min-width: 4;
        margin-right: 1;
    }
    FilePickerModal #current-dir {
        width: 1fr;
        color: $text-muted;
    }
    FilePickerModal FilteredTree {
        height: 1fr;
        border: solid $panel;
    }
    FilePickerModal #selected-path {
        margin-top: 1;
    }
    FilePickerModal #btn-row {
        height: 1;
        margin-top: 1;
        align: right middle;
    }
    FilePickerModal #btn-row Button {
        margin-left: 1;
    }
    """
    BINDINGS = [("escape", "dismiss(None)", "Cancel")]

    def __init__(self, start_path: str = ".", suffixes: frozenset[str] = frozenset({".su2"})) -> None:
        super().__init__()
        self._suffixes = suffixes
        self._root = Path(start_path).parent if start_path else Path.home()
        if not self._root.is_dir():
            self._root = Path.home()
        self._root = self._root.absolute()
        self._selected: str | None = None

    def _valid(self, path: str) -> bool:
        return Path(path).suffix.lower() in self._suffixes

    def compose(self) -> ComposeResult:
        with Vertical(id="dialog"):
            with Horizontal(id="nav-row"):
                yield Button("⮭", id="btn-up", compact=True)
                yield Label(str(self._root), id="current-dir")
            yield FilteredTree(str(self._root), self._suffixes, id="filetree")
            yield Input(placeholder="selected path", id="selected-path", compact=True)
            with Horizontal(id="btn-row"):
                yield Button("Cancel", id="btn-cancel", variant="default", compact=True)
                yield Button("Open", id="btn-open", variant="primary", disabled=True, compact=True)

    @on(Button.Pressed, "#btn-up")
    async def _go_up(self) -> None:
        parent = self._root.parent
        if parent == self._root:
            return
        self._root = parent
        self.query_one("#current-dir", Label).update(str(self._root))
        old_tree = self.query_one("#filetree", DirectoryTree)
        selected_input = self.query_one("#selected-path", Input)
        await old_tree.remove()
        await self.query_one("#dialog").mount(
            FilteredTree(str(self._root), self._suffixes, id="filetree"),
            before=selected_input,
        )

    @on(DirectoryTree.FileSelected)
    def _on_file_selected(self, event: DirectoryTree.FileSelected) -> None:
        path = str(event.path)
        self._selected = path
        self.query_one("#selected-path", Input).value = path
        self.query_one("#btn-open", Button).disabled = not self._valid(path)

    @on(Input.Changed, "#selected-path")
    def _on_path_typed(self, event: Input.Changed) -> None:
        path = event.value.strip()
        self._selected = path if path else None
        self.query_one("#btn-open", Button).disabled = not self._valid(path)

    @on(Input.Submitted, "#selected-path")
    def _on_path_submitted(self, _: Input.Submitted) -> None:
        if self._selected and self._valid(self._selected):
            self.dismiss(self._selected)

    @on(Button.Pressed, "#btn-open")
    def _open(self) -> None:
        self.dismiss(self._selected)

    @on(Button.Pressed, "#btn-cancel")
    def _cancel(self) -> None:
        self.dismiss(None)
