# nonloccfg

Терминальный редактор конфигурационных файлов для конечно-элементного солвера [NonLocFEM](https://github.com/VIXER13/nonlocal).
Позволяет создавать и редактировать JSON-конфиги через удобный TUI вместо ручной правки JSON.

![Python 3.12+](https://img.shields.io/badge/python-3.12%2B-blue)

## Возможности

- **Вкладка Task** — размерность задачи, тип (thermal / mechanical / thermomechanical), зависимость от времени, путь к сетке, граничные условия, материалы
- **Вкладка Materials DB** — локальная база данных материалов с нечётким поиском; заполняет поля из ссылки и блокирует их от случайного изменения
- **Вкладка Options** — пути сохранения результатов, путь к выходному конфигу
- Импорт существующих конфигов по **Ctrl+O**, сохранение по **Ctrl+S**

## Требования

- Python 3.12+
- [uv](https://docs.astral.sh/uv/) (рекомендуется) или pip

## Установка

```bash
uv sync
```

## Запуск

```bash
uv run python main.py
```

| Клавиша | Действие |
|---------|----------|
| `Ctrl+S` | Сохранить конфиг по пути из Options → config output |
| `Ctrl+O` | Открыть и загрузить существующий конфиг |
| `Q` | Выйти |

## Формат конфига

```json
{
    "task": {
        "dimension": 2,
        "problem": "thermomechanical",
        "time_dependency": false
    },
    "mesh": { "path": "mesh/rect.su2" },
    "thermal_boundaries": {
        "Wall": { "kind": "temperature", "temperature": 300.0 }
    },
    "mechanical_boundaries": {
        "Support": [{ "displacement": 0.0 }, { "displacement": 0.0 }]
    },
    "materials": {
        "Body": {
            "physical": {
                "conductivity": 50.0,
                "density": 7800.0,
                "youngs_modulus": 200e9,
                "poissons_ratio": 0.3
            },
            "thermal_model":    { "local_weight": 0.5, "nonlocal_radius": 0.1 },
            "mechanical_model": { "local_weight": 0.5, "nonlocal_radius": 0.1 }
        }
    },
    "save": { "folder": "./results", "vtk": "solution" }
}
```

Для задач типа thermal или mechanical модель материала записывается под ключом `"model"` вместо `"thermal_model"` / `"mechanical_model"`.

## База материалов

База хранится в файле `materials_db.json` в корне проекта. При первом запуске заполняется тремя материалами: Steel, Aluminum, Copper.

**Флоу работы со ссылкой на материал из базы:**

1. Откройте вкладку **Materials DB** и создайте или отредактируйте материал.
2. На вкладке **Task**, после загрузки сетки, нажмите **DB** рядом с нужным слотом материала — поля заполнятся автоматически и заблокируются.
3. Чтобы изменить значения: нажмите **unlink** для отвязки ссылки (значения при этом сохраняются), затем редактируйте.
4. При сохранении (`Ctrl+S`) имя ссылки записывается в конфиг как `"$db_ref"` и восстанавливается при следующем импорте.

## Структура проекта

```
nonloccfg/
├── src/nonloccfg/
│   ├── app.py              # NonlocCfgApp, точка входа
│   ├── config_helpers.py   # _parse_value / _format_value / вспомогательные функции
│   ├── db.py               # синглтон MaterialsDB
│   ├── materials_db.py     # класс MaterialsDB и нечёткий поиск
│   ├── messages.py         # Textual-сообщения (ProblemChanged, MeshLoaded, TagReleased)
│   ├── su2_parser.py       # парсер маркеров SU2-сетки
│   └── widgets/
│       ├── boundaries.py        # ThermalBCWidget, MechanicalBCWidget, BoundariesSection
│       ├── file_picker.py       # FilePickerModal
│       ├── material_picker.py   # MaterialPickerModal (нечёткий поиск + превью)
│       ├── materials.py         # ModelFields, MaterialWidget, MaterialsSection
│       ├── materials_db_tab.py  # MaterialsDBTab
│       ├── options_tab.py       # SaveSection, AuxiliarySection, OutputSection, OptionsTab
│       └── task_tab.py          # TaskSection, MeshSection, TimeSection, TaskTab
├── tests/
│   ├── test_config_helpers.py   # юнит-тесты _parse_value / _format_value
│   ├── test_materials_db.py     # юнит-тесты MaterialsDB и нечёткого поиска
│   ├── test_messages.py         # smoke-тесты сообщений
│   ├── test_su2_parser.py       # юнит-тесты парсера SU2
│   ├── test_widgets.py          # round-trip тесты collect/load виджетов
│   └── test_integration.py      # интеграционные тесты через Textual Pilot
├── main.py                 # тонкая точка входа
├── materials_db.json       # локальная база материалов
└── pyproject.toml
```

## Тесты

```bash
uv run pytest
```
