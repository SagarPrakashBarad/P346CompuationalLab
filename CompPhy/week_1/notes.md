## 4th January 2024

```
Classnotes
```

### Solution of non-linear equations

**one variable equation**   

$cos (x) - x^3 = 0$ 
This can be solved using bisection, regular falsi and newton methods.

**multi-dimensional equation**
 $xe^y -1 = 0 \\ y^2 - x^2 -1 = 0$
Newton raphson and fixed point analysis.


### Numerical integration

$I =  \int_a^b f(x) dx$ --> $I_n$

Considering the $f(x)$ is smooth and have no sigularity in a range then you can use $f(x) = a_0 + a_1 x + a_2 x^2$.

On expanding your integration on legedre polynomials you can the gaussian quadrature.

`plot goes here`

#### Steps for integration

Reformulate the integration as a summation;

$ I = \int_a^b f(x) dx = \sum_{n=1}^N w(x_n) f(x_n) $

* $w(x_n) : weight function$
* $\frac{b-a}{h} = N$

`check that is it laguere or legendre polynomial`


$ I = \int_a^b f(x) dx  = \int_a^b p(x) dx$

The idea is we can approximate f(x) as a polynomial but for gaussian quadrature we need to expand as weighted function of legendre polynomial.


$ I = \int_a^b  [f(x_0)l_0(x) + f(x_1)l_1(x) + f(x_1)l_1(x) + f(x_{N-1})l_{N-1}(x)] dx$


`plot goes here`


