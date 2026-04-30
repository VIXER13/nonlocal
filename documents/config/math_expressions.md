# Интерпретатор математических выражений

В программе реализован интерпретатор математических выражений, который принимает выражения в виде строки следующего вида

> "x y : x^2 + sin(y) + 5.0"

Строка с выражением имеет разделительный символ **:**, слева от символа перечислены используемые в выражении переменные, а справа интерпретируемое выражение. Порядок переменных определяет их индексацию внутри программы из-за чего важно соблюдать контракт, который предписывает [структура конфигурационного файла](config_structure.md) для задания тех или иных величин. Количество переменных и их имена не имеют ограничений, за исключением использования зарезервированных функций, операторов и разделительных символом. Все вычисления происходят в арифметике с плавающей точкой, в том числе все числа заданные литералами интерпретируются как числа с плавающей точкой.

## Основной функционал

Большинство функций повторяют стандартную библиотеку математических функций [cmath](https://en.cppreference.com/w/cpp/header/cmath.html) (такие функции имеют гиперссылки).

### Числовые литералы

- Целочисленная запись (при вычислениях будет сконвертировано в число с плавающей точкой);
- Запись с плавающей точкой (например: 0.5, 10.4, 63., .04);
- Экспоненциальная форма записи (например: 0.1e2, 1.2e-3, 3E+4, 5e0).

### Унарные операторы
- Унарный минус &mdash; **~**.

### Бинарные операторы
- Плюс &mdash; **+**;
- Минус &mdash; **-**;
- Умножение &mdash; **\***;
- Деление &mdash; **/**;
- Возведение в степень &mdash; **^**.

### Функции одного переменного
- Функция знака (сигнум) &mdash; **sign**;
- Абсолютное значение &mdash; **[abs](https://en.cppreference.com/cpp/numeric/math/fabs)**;
- Возведение в квадрат &mdash; **sqr**;
- Извлечение квадратного корня &mdash; **[sqrt](https://en.cppreference.com/cpp/numeric/math/sqrt)**;
- Извлечение кубического корня &mdash; **[cbrt](https://en.cppreference.com/cpp/numeric/math/cbrt)**;
- Синус &mdash; **[sin](https://en.cppreference.com/cpp/numeric/math/sin)**;
- Косинус &mdash; **[cos](https://en.cppreference.com/cpp/numeric/math/cos)**;
- Тангенс &mdash; **[tan](https://en.cppreference.com/cpp/numeric/math/tan)**;
- Арксинус &mdash; **[asin](https://en.cppreference.com/cpp/numeric/math/asin)**;
- Арккосинус &mdash; **[acos](https://en.cppreference.com/cpp/numeric/math/acos)**;
- Арктангенс &mdash; **[atan](https://en.cppreference.com/cpp/numeric/math/atan)**;
- Гиперболический синус &mdash; **[sinh](https://en.cppreference.com/cpp/numeric/math/sinh)**;
- Гиперболический косинус &mdash; **[cosh](https://en.cppreference.com/cpp/numeric/math/cosh)**;
- Гиперболический тангенс &mdash; **[tanh](https://en.cppreference.com/cpp/numeric/math/tanh)**;
- Гиперболический арксинус &mdash; **[asinh](https://en.cppreference.com/cpp/numeric/math/asinh)**;
- Гиперболический арккосинусо &mdash; **[acosh](https://en.cppreference.com/cpp/numeric/math/acosh)**;
- Гиперболический арктангенс &mdash; **[atanh](https://en.cppreference.com/cpp/numeric/math/atanh)**;
- Экспонента с основанием e &mdash; **[exp](https://en.cppreference.com/cpp/numeric/math/exp)**;
- Экспонента с основанием 2 &mdash; **[exp2](https://en.cppreference.com/cpp/numeric/math/exp2)**;
- Экспонента с основанием e минус 1 (т.е. $e^x-1$) &mdash; **[expm1](https://en.cppreference.com/cpp/numeric/math/expm1)**;
- Логарифм с основанием e &mdash; **[log](https://en.cppreference.com/cpp/numeric/math/log)**;
- Логарифм с основанием 2 &mdash; **[log2](https://en.cppreference.com/cpp/numeric/math/log2)**;
- Логарифм с основанием 10 &mdash; **[log10](https://en.cppreference.com/cpp/numeric/math/log10)**;
- Логарифм с основанием e (т.е. $\ln(x + 1)$) &mdash; **[log1p](https://en.cppreference.com/cpp/numeric/math/log1p)**;
- Функция ошибок &mdash; **[erf](https://en.cppreference.com/cpp/numeric/math/erf)**;
- Дополнительная функция ошибок &mdash; **[erfc](https://en.cppreference.com/cpp/numeric/math/erfc)**;
- Гамма функция &mdash; **[tgamma](https://en.cppreference.com/cpp/numeric/math/tgamma)**;
- Натуральный логарифм от гамма функции &mdash; **[lgamma](https://en.cppreference.com/cpp/numeric/math/lgamma)**;
- Вычисление ближайшего целого не меньше данного &mdash; **[ceil](https://en.cppreference.com/cpp/numeric/math/ceil)**;
- Вычисление ближайшего целого не больше данного &mdash; **[floor](https://en.cppreference.com/cpp/numeric/math/floor)**;
- Вычисление ближайшего целого не больше данного по абсолютному значению, с сохранением знака &mdash; **[trunc](https://en.cppreference.com/cpp/numeric/math/trunc)**;
- Вычисление целого путём округления по стандартным правилам округления &mdash; **[round](https://en.cppreference.com/cpp/numeric/math/round)**;

### Функции двух переменных
- Возведение в степень &mdash; **[pow](https://en.cppreference.com/cpp/numeric/math/pow)**;
- Вычисление угла по отношению координаты y к координате x &mdash; **[atan2](https://en.cppreference.com/cpp/numeric/math/atan2)**;
- Вычисление длины гипотенузы &mdash; **[hypot](https://en.cppreference.com/cpp/numeric/math/hypot)**;
- Вычисление остатка от деления &mdash; **[fmod](https://en.cppreference.com/cpp/numeric/math/fmod)**;
- Максимум &mdash; **[max](https://en.cppreference.com/cpp/numeric/math/fmax)**;
- Минимум &mdash; **[min](https://en.cppreference.com/cpp/numeric/math/fmin)**.
