# Основные соотношения нелокальной теории

## Линейный интегральный нелокальный оператор
$$
\mathcal{N} [f(\boldsymbol{x})] = 
	p_1 f(\boldsymbol{x}) + 
	p_2 \int\limits_{S'(\boldsymbol{x}') \cap S} 
		\varphi(\boldsymbol{x}, \boldsymbol{x}') f(\boldsymbol{x}')
	dS'(\boldsymbol{x}'),
	\quad
	\boldsymbol{x}' \in S'(\boldsymbol{x}).
$$
- Здесь $f(\boldsymbol{x})$ &mdash; выражение, описывающее сохраняющуюся физическую субстанцию;
- $p_1 > 0$ и $p_2 \geqslant 0$ &mdash; весовые параметры модели такие, что $p_1 + p_2 = 1$;
- $\varphi$ &mdash; функция нелокального влияния, нормированная положительная монотонно убывающая функция в области $S'(\boldsymbol{x})$;
- $\boldsymbol{x}'$ &mdash; точка в области $S'(\boldsymbol{x})$, в которой вычисляется влияние на величины находящиеся в точке $\boldsymbol{x}$;
- $S'(\boldsymbol{x})$ &mdash; область нелокального влияния с центром в точке $\boldsymbol{x} \in S$;
- $S$ &mdash; область занимаемая рассматриваемым телом.

## Функции нелокального влияния

### Метрическая функция
$$
\rho_n(\boldsymbol{x}, \boldsymbol{x}') = 
	\left(
		\left| \dfrac{x_1 - x_1'}{r_1} \right|^n +
		\left| \dfrac{x_2 - x_2'}{r_2} \right|^n
	\right)^{\dfrac{1}{n}}.
$$
При $n \rightarrow \infty$ метрическая функция принимает вид
$$
\rho_{\infty} (\boldsymbol{x}, \boldsymbol{x}') = 
	\max \left( 
		\left| \dfrac{x_1 - x_1'}{r_1} \right|,
		\left| \dfrac{x_2 - x_2'}{r_2} \right|
	\right),
$$
- $r_i > 0$ &mdash; длины главных полуосей;
- $n > 0$ &mdash; метрический параметр.

### Полиномиальное семейство функций нелокального влияния
$$
\varphi_{p,q}^{P}(\boldsymbol{x}, \boldsymbol{x}') =
	\begin{cases}
		A(1 - \rho_n(\boldsymbol{x}, \boldsymbol{x}')^p)^q, \quad &\rho_n(\boldsymbol{x}, \boldsymbol{x}') \leqslant 1, \\
		0, &\rho_n(\boldsymbol{x}, \boldsymbol{x}') > 1,
	\end{cases}
$$
- $p > 0$, $q > 0$ &mdash; параметры плотности распределения;
- $A$ &mdash; нормировочный параметр, который равен

$$
A = \dfrac{np}
	{
		4 r_1 r_2 
		\operatorname{B}\left( \dfrac{1}{n}, \dfrac{1}{n} \right) 
		\operatorname{B}\left( \dfrac{2}{p}, q+1 \right)
	};
$$

при $n \rightarrow \infty$ принимает вид

$$
A = \dfrac{p}{8 r_2 r_2 \operatorname{B}\left( \dfrac{2}{p}, q+1 \right)}
$$

### Экспоненциальное семейство функций нелокального влияния
$$
\varphi_{p,q}^{E} (\boldsymbol{x}, \boldsymbol{x}') =
	A \exp \left(-q\rho_n(\boldsymbol{x}, \boldsymbol{x}')^p \right),
$$
- $p > 0$, $q > 0$ &mdash; параметры плотности распределения;
- $A$ &mdash; нормировочный параметр, который равен

$$
A = 
	\dfrac
	{
		4^{\frac{1}{n}} n p q^{\frac{2}{p}}
	}
	{
		8 r_1 r_2 \operatorname{B} \left( \dfrac{1}{2}, \dfrac{1}{n} \right) \operatorname{\Gamma} \left( \dfrac{2}{p} \right)
	},
$$

при $n \rightarrow \infty$ принимает вид

$$
A = \dfrac{p q^{\frac{2}{p}}}{8 r_1 r_2 \Gamma \left( \dfrac{2}{p} \right)}.
$$

