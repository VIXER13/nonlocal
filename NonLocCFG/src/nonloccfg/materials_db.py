from __future__ import annotations

import json
from pathlib import Path


def _fuzzy_score(query: str, name: str) -> int:
    """Higher = better match. 0 = no match."""
    if not query:
        return 1
    q, n = query.lower(), name.lower()
    if q in n:
        return 100 - n.index(q)
    qi = 0
    for ch in n:
        if qi < len(q) and ch == q[qi]:
            qi += 1
    return max(1, 50 - len(n)) if qi == len(q) else 0


class MaterialsDB:
    FIELDS: tuple[str, ...] = (
        "conductivity", "capacity", "density", "relaxation_time",
        "youngs_modulus", "poissons_ratio", "shear_modulus", "thermal_expansion",
    )
    MODEL_FIELDS: tuple[str, ...] = (
        "local_weight", "nonlocal_radius", "search_radius",
        "distance", "influence", "n", "p", "q",
    )

    def __init__(self, path: Path) -> None:
        self._path = path
        self._data: dict[str, dict] = {}
        self._load()

    def _load(self) -> None:
        if self._path.exists():
            try:
                self._data = json.loads(self._path.read_text(encoding="utf-8"))
            except Exception:
                self._data = {}

    def flush(self) -> None:
        self._path.write_text(
            json.dumps(self._data, indent=4, ensure_ascii=False),
            encoding="utf-8",
        )

    def all_names(self) -> list[str]:
        return sorted(self._data.keys())

    def get(self, name: str) -> dict | None:
        return dict(self._data[name]) if name in self._data else None

    def upsert(self, name: str, data: dict) -> None:
        self._data[name] = {k: v for k, v in data.items() if v is not None}
        self.flush()

    def delete(self, name: str) -> None:
        self._data.pop(name, None)
        self.flush()

    def search(self, query: str) -> list[str]:
        scored = [(n, _fuzzy_score(query, n)) for n in self._data]
        return [n for n, s in sorted(scored, key=lambda x: -x[1]) if s > 0]
