# NonLocal Finite Element Method (NonLocFEM)

**Healthchecks**

[![gcc-13](https://github.com/VIXER13/nonlocal/actions/workflows/gcc-13.yml/badge.svg?branch=dev)](https://github.com/VIXER13/nonlocal/actions/workflows/gcc-13.yml)
[![gcc-14](https://github.com/VIXER13/nonlocal/actions/workflows/gcc-14.yml/badge.svg?branch=dev)](https://github.com/VIXER13/nonlocal/actions/workflows/gcc-14.yml)
[![gcc-15](https://github.com/VIXER13/nonlocal/actions/workflows/gcc-15.yml/badge.svg?branch=dev)](https://github.com/VIXER13/nonlocal/actions/workflows/gcc-15.yml)

[![clang-20](https://github.com/VIXER13/nonlocal/actions/workflows/clang-20.yml/badge.svg?branch=dev)](https://github.com/VIXER13/nonlocal/actions/workflows/clang-20.yml)
[![clang-21](https://github.com/VIXER13/nonlocal/actions/workflows/clang-21.yml/badge.svg?branch=dev)](https://github.com/VIXER13/nonlocal/actions/workflows/clang-21.yml)
[![clang-22](https://github.com/VIXER13/nonlocal/actions/workflows/clang-22.yml/badge.svg?branch=dev)](https://github.com/VIXER13/nonlocal/actions/workflows/clang-22.yml)

Конечно-элементный программный комплекс NonLocFEM предназначен для эффективного решения задач термоупругости с учётом [пространственной нелокальности](./documents/theory/nonlocal_operator.md) на многопроцессорных электронно-вычислительных машинах с общей и распределённой памятью. К основным особенностям программного комплекса стоит отнести возможность решения одномерных и двумерных задач стационарной и нестационарной [теплопроводности](./documents/theory/thermal.md), и задачи [статики](./documents/theory/mechanical.md), для однородных и композитных материалов на неструктурированных сетках с использованием изопараметрических конечных элементов произвольного порядка и формы.
 
В программном комплексе реализована гибкая система параметризации задачи при помощи [конфигурационных файлов](./documents/config/config_structure.md), содержащих структуры в формате JSON и реализован [интерпретатор математических выражений](./documents/config/math_expressions.md). С помощью данных инструментов возможно определить варианты расчёта, задать граничные условия и параметры материалов. В программном комплексе реализована поддержка расчётов с граничными условиями, определяющими температуру, плотность теплового потока, температурное излучение и их комбинации на границах рассматриваемых областей. Так же существует возможность определить свойства среды, такие как изотропность, ортотропность и анизотропность. Коэффициенты тензора теплопроводности могут содержать в себе константы или выражения с зависимостью от пространственных переменных. В одномерных расчётах нестационарной теплопроводности реализована возможность учесть конечную скорость распространения теплового потока.

**Сборка**
Для сборки необходимы необходимы следующие программы и пакеты: gcc/g++ >= 13 или clang/clang++ >= 18, CMake >= 3.16, conan > 2.0 После их установки необходимо последовательно выполнить следующие команды в командной строке
```bash
make build
```

Команда make соберёт проекты NonLocFEM и unit_tests.

Возможно явно указать нужный компилятор (gcc или clang[по умолчанию])
```bash
make build COMPILER=gcc
```

Версию компилятора можно указать через полное имя. Если версия не задана, используется версия и исполняемые файлы из профиля Conan:
```bash
make build COMPILER=gcc-15
```

Возможно явно указать количество потоков (используются все доступные[по умолчанию])
```bash
make build THREADS=4
```

**Запуск**
Перед первым запуском настоятельно рекомендуется запустить юнит тесты, чтобы убедиться, что основные компоненты программы работают корректно
```bash
make run-tests
```

Для запуска программы NonLocFEM необходимо указать путь к исполняемому файлу и в качестве аргумента передать путь к .json файлу с конфигурацией запуска, например
```bash
./build/NonLocFEM/NonLocFEM ./documents/config/examples/thermal_stationary_1d.json
```

Дополнительно можно указать уровень логирования с помощью аргумента `--log-level`:
```bash
./build/NonLocFEM/NonLocFEM ./documents/config/examples/thermal_stationary_1d.json --log-level=Debug
```

Доступные уровни логирования (в порядке возрастания детализации):
- `off` - отключение логирования
- `error` - ошибки
- `warning` - предупреждения
- `info` - информационные сообщения (по умолчанию)
- `debug` - отладочная информация
- `trace` - трассировка выполнения

Уровень логирования регистронезависимый, например `--log-level=DEBUG` или `--log-level=Info` будут работать корректно.

**Разработка**
- Рекомендуется держать файлы CMakeLists.txt в форматированном состоянии. Для этого удобно использовать gersemi
  ```
  pip3 install gersemi
  find . -type f -name CMakeLists.txt -exec python3 -m  gersemi -i --definitions CMakeTools/EmbeddedFiles.cmake {} +
  ```

## SAST Tools

[PVS-Studio](https://pvs-studio.com/en/pvs-studio/?utm_source=website&utm_medium=github&utm_campaign=open_source) - static code analyzer for Enterprise (C, C++, C#, Go, and Java) and Web (JS and TS) development.
