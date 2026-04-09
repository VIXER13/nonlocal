# Теплопроводность

## Уравнение

Нестационарное уравнение теплопроводности

$$c_V \dfrac{\partial T}{\partial t} = \nabla \cdot \boldsymbol{q} - q_V,$$

cтационарное уравнение теплопроводности

$$\nabla \cdot \boldsymbol{q} = q_V,$$

- $T = T(\boldsymbol{x})$ &mdash; поле температуры;
- $t$ &mdash; время;
- $\nabla = \partial/\partial x_i \boldsymbol{e}_i$ &mdash; дифференциальный оператор набла;
- $\boldsymbol{x} = x_i \boldsymbol{e}_i$ &mdash; вектор пространственной переменной;
- $\boldsymbol{e}_i$ &mdash; единичный орт;
- $c_V$ &mdash; удельная объёмная теплоёмкость;
- $q_V$ &mdash; объёмная плотность мощности внутренних источников и стоков теплоты;
- $\boldsymbol{q} = q_i \boldsymbol{e}_i$ &mdash; вектор плотности теплового потока, который определён через обобщённую гипотезу Био &mdash; Фурье

$$\boldsymbol{q}(\boldsymbol{x}) = \mathcal{N} \left( -\widehat{\boldsymbol{\lambda}} \cdot \nabla T \right),$$

- $\mathcal{N}$ &mdash; [линейный интегральный нелокальный оператор](./nonlocal_operator.md#линейный-интегральный-нелокальный-оператор);
- $\widehat{\boldsymbol{\lambda}} = \lambda_{ij} (\boldsymbol{x}) \boldsymbol{e}_i \boldsymbol{e}_j$ &mdash; [тензор коэффициентов теплопроводности](#тензор-теплопроводности).

## Граничные условия

Рассмотрены 4 варианта граничных условий:
1. $T|_{\Gamma_1} = T _{\Gamma}(\boldsymbol{x})$ &mdash; температурное граничное условие (первого рода);
2. $\boldsymbol{n} \cdot \boldsymbol{q}|_{\Gamma_2} = f(\boldsymbol{x})$ &mdash; потоковое граничное условие (второго рода);
3. $\boldsymbol{n} \cdot \boldsymbol{q}|_{\Gamma_3} = \alpha (T_a(\boldsymbol{x}) - T(\boldsymbol{x}))$ &mdash; условие конвективного теплообмена (третьего рода);
4. $\boldsymbol{n} \cdot \boldsymbol{q}|_{\Gamma_4} = \varepsilon_r \sigma T^4$ &mdash; условие собственного излучения.

Здесь $\bigcup\limits_i \Gamma_i = \partial S$, $\bigcap\limits_i\Gamma_i = \varnothing$, где 
- $\Gamma_i$ &mdash; подобласть границы тела $\partial S$;
- $T_\Gamma$ и $f$ &mdash; функции определяющие температуру и величину теплового потока на границе $\Gamma_1$ и $\Gamma_2$ соответственно;
- $\boldsymbol{n}$ &mdash; вектор внешней нормали;
- $\alpha$ &mdash; коэффициент конвективного теплообмена с внешней средой;
- $T_a$ &mdash; температура внешней среды вблизи границы $\Gamma_3$;
- $\varepsilon_r$ &mdash; коэффициент излучения;
- $\sigma = 5.67036713 \cdot 10^{-8} \text{Вт} \cdot \text{м}^{-2} \cdot \text{К}^{-4}$ &mdash; постоянная Стефана &mdash; Больцмана.

Для моделирования идеального теплового контакта используются условия следующего вида

$$
T_i|_{\Gamma_{ij}} = T_j|_{\Gamma_{ij}},
\quad
\boldsymbol{n} \cdot \boldsymbol{q}_i|_{\Gamma_{ij}} = \boldsymbol{n} \cdot \boldsymbol{q}_j|_{\Gamma_{ij}},
$$

где $\Gamma_{ij}$ &mdash; граница контакта материалов под номерами $i$ и $j$.

## Тензор теплопроводности

### Одномерный случай

В одномерном случае тензор теплопроводности $\widehat{\boldsymbol{\lambda}}$ представлен в виде скалярной величины $\lambda$.

### Двумерный случай

#### Изотропный материал

В изотропном случае тензор теплопроводности $\widehat{\boldsymbol{\lambda}}$ представлен в виде скалярной величины $\lambda$.

#### Ортотропный материал

В ортотропном случае тензор теплопроводности $\widehat{\boldsymbol{\lambda}}$ имеет две независимых компоненты

$$
\widehat{\boldsymbol{\lambda}} = 
\begin{pmatrix}
    \lambda_{11} & 0 \\
    0 & \lambda_{22}
\end{pmatrix}
$$

#### Анизотропный материал

В анизотропном случае тензор теплопроводности $\widehat{\boldsymbol{\lambda}}$ имеет три независимых компоненты

$$
\widehat{\boldsymbol{\lambda}} = 
\begin{pmatrix}
    \lambda_{11} & \lambda_{12} \\
    \lambda_{12} & \lambda_{22}
\end{pmatrix}
$$
