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
- $\widehat{\text{\textbf{C}}} = C_{ijkl} \boldsymbol{e}_i \otimes \boldsymbol{e}_j \otimes \boldsymbol{e}_k \otimes \boldsymbol{e}_l$ &mdash; [тензор коэффициентов упругости](#тензор-коэффициентов-упругости);
- $\widehat{\boldsymbol{\alpha}}^T = \alpha_{ij}^T \boldsymbol{e}_i \otimes \boldsymbol{e}_j$ &mdash; [тензор температурных коэффициентов линейного расширения](#тензор-температурных-коэффициентов-линейного-расширения);
- $\Delta T = T - T_0$ &mdash; разница между текущим распределением температуры $T$ и распределением $T_0$ при котором отсутствуют температурные деформации;
- $\widehat{\varepsilon} = \varepsilon_{ij} \boldsymbol{e}_i \otimes \boldsymbol{e}_j$ &mdash; тензор деформации, считаем, что деформации достаточно малы, поэтому для определения компонент тензора деформации $\widehat{\boldsymbol{\varepsilon}}$ воспользуемся соотношением Коши

$$
\widehat{\boldsymbol{\varepsilon}} = 
	\dfrac{\nabla \boldsymbol{u} + (\nabla \boldsymbol{u})^T}{2} = 
	\dfrac{u_{i, j} + u_{j, i}}{2} \boldsymbol{e}_i \otimes \boldsymbol{e}_j,
$$

- $\boldsymbol{u} = u_i \boldsymbol{e}_i$ &mdash; вектор перемещения.

## Граничные условия

Рассмотрены граничные условия двух видов:
1. $\boldsymbol{u}|_{\Gamma_1} = \boldsymbol{d} (\boldsymbol{x})$ &mdash; кинематические граничные условия (первого рода);
2. $\boldsymbol{n} \cdot \widehat{\boldsymbol{\sigma}}|_{\Gamma_2} = \boldsymbol{p} (\boldsymbol{x})$ &mdash; силовые граничные условия (второго рода);

**Примечание**: Так же возможны смешанные условия, когда по одной оси задано кинематическое условие, а по второй силовое.

Здесь $\Gamma_1 \cup \Gamma_2 = \partial S$, $\Gamma_1 \cap \Gamma_2 = \varnothing$, где
- $\Gamma_i$ &mdash; подобласть границы тела $\partial S$;
- $\boldsymbol{d} (\boldsymbol{x}) = d_i (\boldsymbol{x}) \boldsymbol{e}_i$ &mdash; вектор перемещений на границе $\Gamma_1$;
- $\boldsymbol{p} (\boldsymbol{x}) = p_i (\boldsymbol{x}) \boldsymbol{e}_i$ &mdash; вектор плотности поверхностностного нагружения на границе $\Gamma_2$;
- $\boldsymbol{n}$ &mdash; вектор внешней нормали.

Для моделирования идеального механического контакта используются условия следующего вида

$$
\boldsymbol{u}_i|_{\Gamma_{ij}} = \boldsymbol{u}_j|_{\Gamma_{ij}},
\quad
\boldsymbol{n} \cdot \widehat{\boldsymbol{\sigma}}_i|_{\Gamma_{ij}} = \boldsymbol{n} \cdot \widehat{\boldsymbol{\sigma}}_j|_{\Gamma_{ij}},
$$

где $\Gamma_{ij}$ &mdash; граница контакта материалов под номерами $i$ и $j$.

## Тензор коэффициентов упругости

### Двумерный случай

#### Изотропный материал

В изотропном случае параметры материала могут быть заданы определением двух упругих переменных:
- $E$ &mdash; модуль Юнга;
- $\nu$ &mdash; коэффициент Пуассона.

Матрица Гука имеет вид

$$
D = 
\begin{pmatrix}
	d_{11} & d_{12} & 0 \\
	d_{12} & d_{11} & 0 \\
	     0 &      0 & d_{66}
\end{pmatrix} = \dfrac{E}{1 - \nu^2}
\begin{pmatrix}
	  1 & \nu & 0 \\
	\nu &   1 & 0 \\
	  0 &   0 & \dfrac{1 - \nu}{2}
\end{pmatrix}
$$

#### Ортотропный материал

В ортотропном случае параметры материала могут быть определены пятью упругими переменными:
- $E_1$, $E_2$ &mdash; модули Юнга вдоль главных направлений;
- $\nu_{12}$, $\nu_{21}$ &mdash; коэффициенты Пуассона;
- $G$ &mdash; модуль сдвига.

Матрица Гука имеет вид

$$
D = 
\begin{pmatrix}
	d_{11} & d_{12} & 0 \\
	d_{12} & d_{22} & 0 \\
	     0 &      0 & d_{66}
\end{pmatrix} = \dfrac{1}{1 - \nu_{12} \nu_{21}}
\begin{pmatrix}
	         E_1 & \nu_{12} E_1 & 0 \\
	\nu_{21} E_2 &          E_2 & 0 \\
	           0 &            0 & G (1 - \nu_{12} \nu_{21})
\end{pmatrix}
$$

**Примечание**: учитывая симметрию матрицы Гука ($\nu_{12} E_1 = \nu_{21} E_2$), таким образом определение количество необходимых параметров материала можно сократить до четырёх.

#### Анизотропный материал

В анизотропном случае параметры материала могут быть определены аналогично [ортотропному](#ортотропный-материал) случаю c дополнительным параметром поворота $\theta$.

Матрица Гука имеет вид

$$
D = 
\begin{pmatrix}
	d_{11} & d_{12} & d_{16} \\
	d_{12} & d_{22} & d_{26} \\
	d_{16} & d_{26} & d_{66}
\end{pmatrix} = T_1^{-1}(\theta)
\begin{pmatrix}
	d_{11} & d_{12} & 0 \\
	d_{12} & d_{22} & 0 \\
	     0 &      0 & d_{66}
\end{pmatrix}
T_2(\theta),
$$

где матрицы поворота определены следующим образом

$$
T_1 (\theta) =
\begin{pmatrix}
	\cos^2 \theta & \sin^2 \theta &  2 \sin \theta \cos \theta \\
	\sin^2 \theta & \cos^2 \theta & -2 \sin \theta \cos \theta \\
	-\sin \theta \cos \theta & \sin \theta \cos \theta & \cos^2 \theta - \sin^2 \theta
\end{pmatrix}
$$

$$
T_2 (\theta) =
\begin{pmatrix}
	\cos^2 \theta & \sin^2 \theta &  \sin \theta \cos \theta \\
	\sin^2 \theta & \cos^2 \theta & -\sin \theta \cos \theta \\
	-2 \sin \theta \cos \theta & 2 \sin \theta \cos \theta & \cos^2 \theta - \sin^2 \theta
\end{pmatrix}
$$

## Тензор температурных коэффициентов линейного расширения

### Двумерный случай

#### Изотропный материал

В изотропном случае тензор теплопроводности $\widehat{\boldsymbol{\alpha}}^T$ представлен в виде скалярной величины $\alpha^T$.

#### Ортотропный материал

В ортотропном случае тензор теплопроводности $\widehat{\boldsymbol{\alpha}}^T$ имеет две независимых компоненты

$$
\widehat{\boldsymbol{\alpha}}^T = 
\begin{pmatrix}
    \alpha_{11}^T & 0 \\
    0 & \alpha_{22}^T
\end{pmatrix}
$$

#### Анизотропный материал

В анизотропном случае тензор теплопроводности $\widehat{\boldsymbol{\alpha}}^T$ имеет три независимых компоненты

$$
\widehat{\boldsymbol{\alpha}}^T = 
\begin{pmatrix}
    \alpha_{11}^T & \alpha_{12}^T \\
    \alpha_{12}^T & \alpha_{22}^T
\end{pmatrix}
$$
