from __future__ import annotations

import json

from textual.widgets import Input, RadioSet, RadioButton, Select


def _format_value(v) -> str:
    """Convert a Python value back to an Input string."""
    if isinstance(v, list):
        return json.dumps(v)
    return str(v)


def _parse_value(raw: str) -> int | float | list | str | None:
    """Convert a raw input string to int, float, list, or str. Empty string → None."""
    s = raw.strip()
    if not s:
        return None
    if s.startswith("["):
        try:
            parsed = json.loads(s)
            if isinstance(parsed, list):
                return parsed
        except (json.JSONDecodeError, ValueError):
            pass
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        pass
    return s


def _input(widget, id_: str) -> str:
    """Get the value of an Input widget by id, returning '' if absent."""
    try:
        return widget.query_one(f"#{id_}", Input).value
    except Exception:
        return ""


def _radio(widget, id_: str) -> str:
    """Get the plain-text label of the selected RadioButton in a RadioSet."""
    try:
        rs = widget.query_one(f"#{id_}", RadioSet)
        # Use _selected index directly — pressed_button reads RadioButton.value
        # which only updates after a reactive cycle, causing stale reads after
        # programmatic _set_radio calls.
        if rs._selected is None:
            return ""
        buttons = list(rs.query(RadioButton))
        return buttons[rs._selected].label.plain
    except Exception:
        return ""


def _set_radio(widget, id_: str, value: str) -> None:
    """Select the RadioButton whose label matches value."""
    rs = widget.query_one(f"#{id_}", RadioSet)
    for i, btn in enumerate(rs.query(RadioButton)):
        if btn.label.plain == value:
            rs._selected = i
            return


def _select(widget, id_: str) -> str:
    """Get the string value of a Select widget, returning '' if NULL."""
    try:
        sel = widget.query_one(f"#{id_}", Select)
        return "" if sel.value is Select.NULL else str(sel.value)
    except Exception:
        return ""
