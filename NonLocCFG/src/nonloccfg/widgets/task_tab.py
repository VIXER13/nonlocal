from __future__ import annotations

from pathlib import Path

from textual import on
from textual.app import ComposeResult
from textual.containers import Horizontal, VerticalScroll
from textual.widget import Widget
from textual.widgets import Button, Input, Label, RadioButton, RadioSet, Select

from nonloccfg.config_helpers import _format_value, _input, _parse_value, _radio, _select, _set_radio
from nonloccfg.messages import MeshLoaded, ProblemChanged
from nonloccfg.su2_parser import parse_su2_markers
from nonloccfg.widgets.file_picker import FilePickerModal


class TaskSection(Widget):
    DEFAULT_CSS = """
    TaskSection { height: auto; margin-bottom: 1; }
    TaskSection .row { height: 1; align: left middle; padding-left: 2; margin-top: 1; }
    TaskSection .key { color: $text-muted; width: 22; }
    TaskSection .row Input, TaskSection .row Select { width: 1fr; }
    """

    def compose(self) -> ComposeResult:
        yield Label("[bold][#5f87ff]task[/][/]:")
        with Horizontal(classes="row"):
            yield Label("dimension:", classes="key")
            with RadioSet(id="task-dimension"):
                yield RadioButton("1", disabled=True)
                yield RadioButton("2", value=True)
        with Horizontal(classes="row"):
            yield Label("problem:", classes="key")
            yield Select(
                [("thermal", "thermal"), ("mechanical", "mechanical"), ("thermomechanical", "thermomechanical")],
                value="thermal", id="task-problem", compact=True,
            )
        with Horizontal(classes="row"):
            yield Label("time_dependency:", classes="key")
            with RadioSet(id="task-time-dependency"):
                yield RadioButton("false", value=True)
                yield RadioButton("true")

    @on(Select.Changed, "#task-problem")
    def _on_problem(self, event: Select.Changed) -> None:
        if event.value is not Select.NULL:
            self.post_message(ProblemChanged(str(event.value)))

    def collect(self) -> dict:
        return {
            "dimension": int(_radio(self, "task-dimension") or "2"),
            "problem": _select(self, "task-problem") or "thermal",
            "time_dependency": _radio(self, "task-time-dependency") == "true",
        }

    def load(self, data: dict) -> None:
        if "problem" in data:
            self.query_one("#task-problem", Select).value = data["problem"]
        if "time_dependency" in data:
            _set_radio(self, "task-time-dependency", "true" if data["time_dependency"] else "false")


class MeshSection(Widget):
    DEFAULT_CSS = """
    MeshSection { height: auto; margin-bottom: 1; }
    MeshSection .row { height: 1; align: left middle; padding-left: 2; margin-top: 1; }
    MeshSection .key { color: $text-muted; width: 22; }
    MeshSection .row Input { width: 1fr; }
    MeshSection .row Button { width: 10; min-width: 10; }
    MeshSection #mesh-status { color: $text-muted; padding-left: 2; }
    """

    def __init__(self) -> None:
        super().__init__()
        self._current_path: str = "."

    def compose(self) -> ComposeResult:
        yield Label("[bold][#5f87ff]mesh[/][/]:")
        with Horizontal(classes="row"):
            yield Label("path:", classes="key")
            yield Input(placeholder="/path/to/mesh.su2", id="mesh-path-input", compact=True)
            yield Button("Browse", id="mesh-browse", compact=True)
        yield Label("", id="mesh-status")

    @on(Button.Pressed, "#mesh-browse")
    def _browse(self) -> None:
        self.app.push_screen(FilePickerModal(self._current_path), self._on_picked)

    def _on_picked(self, path: str | None) -> None:
        if not path:
            return
        self._current_path = path
        self.query_one("#mesh-path-input", Input).value = path
        self._load(path)

    @on(Input.Submitted, "#mesh-path-input")
    @on(Input.Blurred, "#mesh-path-input")
    def _on_typed(self, _) -> None:
        path = self.query_one("#mesh-path-input", Input).value.strip()
        if path and path != self._current_path:
            self._current_path = path
            self._load(path)

    def _get_dimension(self) -> int:
        rs = self.app.query_one("#task-dimension", RadioSet)
        return int(rs.pressed_button.label.plain)

    def collect(self) -> dict:
        return {"path": _input(self, "mesh-path-input")}

    def load(self, data: dict) -> None:
        path = data.get("path", "")
        if path:
            self._current_path = path
            self.query_one("#mesh-path-input", Input).value = path
            self._load(path)

    def _load(self, path: str) -> None:
        status = self.query_one("#mesh-status", Label)
        dimension = self._get_dimension()
        btags, mtags = parse_su2_markers(path, dimension)
        if not btags and not mtags:
            status.update("  [red]no markers found[/]")
        else:
            status.update(
                f"  [green]+[/] {len(btags) + len(mtags)} markers: "
                f"[dim]{len(btags)} boundaries, {len(mtags)} materials[/]"
            )
            self.post_message(MeshLoaded(btags, mtags))


class TimeSection(Widget):
    DEFAULT_CSS = """
    TimeSection { height: auto; margin-bottom: 1; }
    TimeSection .row { height: 1; align: left middle; padding-left: 2; margin-top: 1; }
    TimeSection .key { color: $text-muted; width: 22; }
    TimeSection .row Input { width: 1fr; }
    """

    def compose(self) -> ComposeResult:
        yield Label("[bold][#5f87ff]time[/][/]:")
        with Horizontal(classes="row"):
            yield Label("time_step:", classes="key")
            yield Input(placeholder="e.g. 0.01", id="time-step", compact=True)
        with Horizontal(classes="row"):
            yield Label("initial_time:", classes="key")
            yield Input(placeholder="0.0  # default", id="time-initial", compact=True)
        with Horizontal(classes="row"):
            yield Label("steps_count:", classes="key")
            yield Input(placeholder="e.g. 100", id="time-steps-count", compact=True)
        with Horizontal(classes="row"):
            yield Label("save_frequency:", classes="key")
            yield Input(placeholder="1  # default", id="time-save-freq", compact=True)

    def collect(self) -> dict:
        d: dict = {}
        for key, wid in (
            ("time_step",      "time-step"),
            ("initial_time",   "time-initial"),
            ("steps_count",    "time-steps-count"),
            ("save_frequency", "time-save-freq"),
        ):
            v = _parse_value(_input(self, wid))
            if v is not None:
                d[key] = v
        return d

    def load(self, data: dict) -> None:
        for key, wid in (
            ("time_step",      "time-step"),
            ("initial_time",   "time-initial"),
            ("steps_count",    "time-steps-count"),
            ("save_frequency", "time-save-freq"),
        ):
            if key in data:
                self.query_one(f"#{wid}", Input).value = _format_value(data[key])


class TaskTab(Widget):
    DEFAULT_CSS = """
    TaskTab { height: 1fr; }
    TaskTab VerticalScroll { padding: 0 2; }
    """

    def compose(self) -> ComposeResult:
        from nonloccfg.widgets.boundaries import BoundariesSection
        from nonloccfg.widgets.materials import MaterialsSection
        with VerticalScroll():
            yield TaskSection()
            yield MeshSection()
            yield TimeSection()
            yield BoundariesSection()
            yield MaterialsSection()

    def on_mount(self) -> None:
        self.query_one(TimeSection).display = False

    @on(RadioSet.Changed, "#task-time-dependency")
    def _toggle_time(self, event: RadioSet.Changed) -> None:
        self.query_one(TimeSection).display = event.pressed.label.plain == "true"

    def on_problem_changed(self, event: ProblemChanged) -> None:
        from nonloccfg.widgets.boundaries import BoundariesSection
        from nonloccfg.widgets.materials import MaterialsSection
        self.query_one(BoundariesSection).show_problem(event.problem)
        self.query_one(MaterialsSection).show_problem(event.problem)

    def on_mesh_loaded(self, event: MeshLoaded) -> None:
        from nonloccfg.widgets.boundaries import BoundariesSection
        from nonloccfg.widgets.materials import MaterialsSection
        self.query_one(BoundariesSection).set_boundary_tags(event.boundary_tags)
        self.query_one(MaterialsSection).set_material_tags(event.material_tags)

    async def load_config(self, cfg: dict, config_dir: Path) -> None:
        from nonloccfg.widgets.boundaries import BoundariesSection
        from nonloccfg.widgets.materials import MaterialsSection
        task = cfg.get("task", {})
        self.query_one(TaskSection).load(task)
        problem = task.get("problem", "thermal")
        time_dep = task.get("time_dependency", False)

        self.query_one(BoundariesSection).show_problem(problem)
        self.query_one(MaterialsSection).show_problem(problem)
        self.query_one(TimeSection).display = time_dep

        mesh_path_raw = cfg.get("mesh", {}).get("path", "")
        if mesh_path_raw:
            mesh_path_abs = (config_dir / mesh_path_raw).resolve()
            if not mesh_path_abs.exists():
                for parent in config_dir.parents:
                    candidate = (parent / mesh_path_raw).resolve()
                    if candidate.exists():
                        mesh_path_abs = candidate
                        break
            mesh_cfg = dict(cfg["mesh"])
            mesh_cfg["path"] = str(mesh_path_abs)
            self.query_one(MeshSection).load(mesh_cfg)
            dimension = task.get("dimension", 2)
            btags, mtags = parse_su2_markers(str(mesh_path_abs), dimension)
            self.query_one(BoundariesSection).set_boundary_tags(btags)
            self.query_one(MaterialsSection).set_material_tags(mtags)

        if time_dep and "time" in cfg:
            self.query_one(TimeSection).load(cfg["time"])

        await self.query_one(BoundariesSection).load(cfg, problem)
        await self.query_one(MaterialsSection).load(cfg.get("materials", {}), problem)

    def build_config(self) -> dict:
        from nonloccfg.widgets.boundaries import BoundariesSection
        from nonloccfg.widgets.materials import MaterialsSection
        task = self.query_one(TaskSection).collect()
        problem = task["problem"]
        cfg: dict = {"task": task}
        cfg["mesh"] = self.query_one(MeshSection).collect()
        if task["time_dependency"]:
            cfg["time"] = self.query_one(TimeSection).collect()
        cfg.update(self.query_one(BoundariesSection).collect(problem))
        cfg["materials"] = self.query_one(MaterialsSection).collect(problem)
        return cfg
