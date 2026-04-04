# Механика

Уравнение равновесия

$$\nabla \cdot \widehat{\boldsymbol{\sigma}} = -\boldsymbol{b},$$

- $\nabla = \partial/\partial x_i \boldsymbol{e}_i$ &mdash; дифференциальный оператор набла;
- $\boldsymbol{x} = x_i \boldsymbol{e}_i$ &mdash; вектор пространственной переменной;
- $\boldsymbol{e}_i$ &mdash; единичный орт;
- $\boldsymbol{b} = b_i \boldsymbol{e}_i$ &mdash; вектор плотности объёмных сил;
- $\widehat{\boldsymbol{\sigma}} = \sigma_{ij} \boldsymbol{e}_i \otimes \boldsymbol{e}_j$ &mdash; тензор напряжений, который определён через обобщённый закон Дюамеля &mdash; Неймана

$$
\widehat{\boldsymbol{\sigma}}(\boldsymbol{x}) =
	\mathcal{N} \left(
		\widehat{\text{\textbf{C}}} \cdot \cdot 
		\left( \widehat{\boldsymbol{\varepsilon}} - \widehat{\boldsymbol{\alpha}}^T \Delta T \right)
	\right),
$$

- $\mathcal{N}$ &mdash; [линейный интегральный нелокальный оператор](./nonlocal_operator.md#линейный-интегральный-нелокальный-оператор);
- $\widehat{\text{\textbf{C}}} = C_{ijkl} \boldsymbol{e}_i \otimes \boldsymbol{e}_j \otimes \boldsymbol{e}_k \otimes \boldsymbol{e}_l$ &mdash; тензор коэффициентов упругости;
- $\widehat{\boldsymbol{\alpha}}^T = \alpha_{ij}^T \boldsymbol{e}_i \otimes \boldsymbol{e}_j$ &mdash; тензор температурных коэффициентов линейного расширения;
- $\Delta T = T - T_0$ &mdash; разница между текущим распределением температуры $T$ и распределением $T_0$ при котором отсутствуют температурные деформации;
- $\widehat{\varepsilon} = \varepsilon_{ij} \boldsymbol{e}_i \otimes \boldsymbol{e}_j$ &mdash; тензор деформации, считаем, что деформации достаточно малы, поэтому для определения компонент тензора деформации $\widehat{\boldsymbol{\varepsilon}}$ воспользуемся соотношением Коши

$$
\widehat{\boldsymbol{\varepsilon}} = 
	\dfrac{\nabla \boldsymbol{u} + (\nabla \boldsymbol{u})^T}{2} = 
	\dfrac{u_{i, j} + u_{j, i}}{2} \boldsymbol{e}_i \otimes \boldsymbol{e}_j,
$$

- $\boldsymbol{u} = u_i \boldsymbol{e}_i$ &mdash; вектор перемещения.
