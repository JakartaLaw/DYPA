# Numerical Integration + Simulation

### Simulation

![Model with continous shocks](assets/markdown-img-paste-20190212133526939.png)

- We find that it's hard to evaluate the model's expectation.
- Taking expectation requires solving af integral. This can be infeasible (f.x. normal distribution don't have a closed form solution for a CDF).
- Solve the problem as a discrete sum:
- ![Discrete sum integral](assets/markdown-img-paste-20190212133917766.png)
- There exists 4 possible solutions:
    - Monte Carlo Integration
    - Equiprobable Integration
    - Gaussian quadrature (Usually this is the best way)
    - Discretization

### Montecarlo integration

- Draw random $S$ numbers $x_i$ from som CDF.
- Calculate: ![Monte Carlo Integration equation](assets/markdown-img-paste-20190212134150500.png)

### Equiprobable points

![Illustration of Equiprobable Points](assets/markdown-img-paste-20190212134334459.png)

Construct equally spaced bins. Do montecarlo integration within each bin and use the weight (ie. The probability of being in given bin), to calculate the Expectation.

![Equiprobable points](assets/markdown-img-paste-20190212134510302.png)

### Gaussian Quadrature

> __There are formulas__ for the sequences of $x_i$ and $\omega_i$ for __exact__ integration of certain polynomials

The formula:

![Gaussian Quadrature](assets/markdown-img-paste-20190212134736278.png)

__Look at algorithm 10 for implementation of model__

### Euler Equation

Assume we not to optimize to functions, where one takes the output of the other: $y=f(x)$ and $z = g(x, y)$. This can be solved with the envelope theorem.

![Envelope theorem](assets/markdown-img-paste-20190212135642591.png)

Using this we can find the optimal relationship between two periods:

![Euler-equation](assets/markdown-img-paste-20190212135724252.png)

This is the same as the Ramsey Model in Macro 3.

### The euler residual

![Euler residual](assets/markdown-img-paste-20190212135839192.png)

A measure of the model (how much the stochastic part obscures the result compared to an un-deterministic model).
