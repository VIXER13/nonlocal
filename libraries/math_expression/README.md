# formula_parser

Contains implementations of run time formula parser

# run_time_parser

Special input format is required 
* pre_infix_notation -> |variables : expression|. Example: |x y z t : x * y * z - t / x + sin(x * y * z)|.
* infix notation (standart) -> x + 5 * (y - z / t), polish notation (prefix) -> x 5 y z t / - * +
  
We convert pre_infix_notation to polish_notation in order to further substituting variables into corresponding formula and its evaluation.

# How to use
## Example №1
```c++
#include "math_expression.h"

int main(int argc, char** argv) {
	auto test = math_expression("v1 v2 v3 p1 p2 p3 m : (p1^2 + p2^2 + p3^2) / (2 * m) + .5 * (v1^2 + v2^2 + v3^2) ");
	assert(test.get_polish_notation() == "p12^p22^+p32^+2m*/.5v12^v22^+v32^+*+" && test.get_n_vars() == 7);
	assert(std::abs(test({ 1, 1, 1, 2, 2, 2, 2 }) - 4.5) < 1e-8);
	
	test = math_expression("x a b : (x - 1)^(a - 1) * (x + 1)^(b + 1)");
	assert(test.get_polish_notation() == "x1-a1-^x1+b1+^*");
	assert(test.get_n_vars() == 3);
	assert(std::abs(test({ 1., 1., 1. }) - 4) < 1e-15);
	assert(std::abs(test({ 1., 100., 100. })) < 1e-15);
	assert(std::abs(test({ -1., 100., 100. })) < 1e-15);
	assert(std::abs(test({ 3., 3., 2. }) - 256) < 1e-15);
	
	test = math_expression("x y t : x * cos(y) / exp(t) + 10");
	assert(test.get_polish_notation() == "xycos*texp/10+");
	assert(test.get_n_vars() == 3);
	assert(std::abs(test({ 2., M_PI, 10. }) - 9.99991) < 1e-5);
	return 0;
}
```

You can use math_expression as ordinary scalar function of vector argument.
