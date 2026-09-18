from __future__ import annotations

from pathlib import Path
from nonloccfg.materials_db import MaterialsDB

# parents[0]=nonloccfg, parents[1]=src, parents[2]=project root
DB_PATH = Path(__file__).parents[2] / "materials_db.json"
db = MaterialsDB(DB_PATH)

# Seed with a few common materials if DB is empty
if not db.all_names():
    db.upsert("Steel", {
        "youngs_modulus": 200e9, "poissons_ratio": 0.3,
        "density": 7800.0, "thermal_expansion": 12e-6,
        "conductivity": 50.0, "capacity": 500.0,
    })
    db.upsert("Aluminum", {
        "youngs_modulus": 70e9, "poissons_ratio": 0.33,
        "density": 2700.0, "thermal_expansion": 23e-6,
        "conductivity": 237.0, "capacity": 900.0,
    })
    db.upsert("Copper", {
        "youngs_modulus": 110e9, "poissons_ratio": 0.34,
        "density": 8960.0, "thermal_expansion": 17e-6,
        "conductivity": 401.0, "capacity": 385.0,
    })
