# Теплопроводность

Нестационарное уравнение теплопроводности

$$c_V \dfrac{\partial T}{\partial t} = \nabla \cdot \boldsymbol{q} - q_V,$$

cтационарное уравнение теплопроводности

$$\nabla \cdot \boldsymbol{q} = q_V,$$

- $T = T(\boldsymbol{x})$ &mdash; поле температуры;
- $t$ &mdash; время;
- $\boldsymbol{x} = x_i \boldsymbol{e}_i$ &mdash; вектор пространственной переменной;
- $\boldsymbol{e}_i$ &mdash; единичный орт;
- $\nabla = \partial/\partial x_i \boldsymbol{e}_i$ &mdash; дифференциальный оператор набла;
- $c_V$ &mdash; удельная объёмная теплоёмкость;
- $q_V$ &mdash; объёмная плотность мощности внутренних источников и стоков теплоты;
- $\boldsymbol{q} = q_i \boldsymbol{e}_i$ &mdash; вектор плотности теплового потока, который определён через обобщённую гипотезу Био &mdash; Фурье

$$\boldsymbol{q}(\boldsymbol{x}) = \mathcal{N} \left( -\widehat{\boldsymbol{\lambda}} \cdot \nabla T \right),$$

- $\mathcal{N}$ &mdash; [линейный интегральный нелокальный оператор](./nonlocal_operator.md#линейный-интегральный-нелокальный-оператор);
- $\widehat{\boldsymbol{\lambda}} = \lambda_{ij} (\boldsymbol{x}) \boldsymbol{e}_i \boldsymbol{e}_j$ тензор коэффициентов теплопроводности.

### Тензор теплопроводности в 1D

В одномерном случае тензор теплопроводности $\widehat{\boldsymbol{\lambda}}$ представлен в виде скалярной величины $\lambda$.

### Тензор теплопроводности в 2D

#### Изотропный случай

В изотропном случае тензор теплопроводности $\widehat{\boldsymbol{\lambda}}$ представлен в виде скалярной величины $\lambda$.

#### Ортотропный случай

В ортотропном случае тензор теплопроводности $\widehat{\boldsymbol{\lambda}}$ имеет две независимых компоненты

$$
\widehat{\boldsymbol{\lambda}} = 
\begin{pmatrix}
    \lambda_{11} & 0 \\
    0 & \lambda_{22}
\end{pmatrix}
$$

#### Анизотропный случай

В ортотропном случае тензор теплопроводности $\widehat{\boldsymbol{\lambda}}$ имеет три независимых компоненты

$$
\widehat{\boldsymbol{\lambda}} = 
\begin{pmatrix}
    \lambda_{11} & \lambda_{12} \\
    \lambda_{12} & \lambda_{22}
\end{pmatrix}
$$
