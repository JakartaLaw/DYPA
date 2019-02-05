# Introduction

![Simplest model](assets/markdown-img-paste-20190203132251674.png)

Mulitple algorithms proposed for solution

### Backwards Induction

Solving the problem by iteratively solving some subset, through recursion:

Illustrated as:
![Backwards induction](assets/markdown-img-paste-20190203133351941.png)

The concept of maximization of $V_{T-2}$, by maximizing $V_{t-1}$ and so on and so forth.

Two important algorithms:

![Find $V$ algorithm](assets/markdown-img-paste-20190203134415209.png)

![FInd all $V$ algorithm](assets/markdown-img-paste-20190203134449286.png)


Important concepts are:

- States (number of goods $M_t$)
- Choices (Consumption $C_t$)
- Payoff function (utility , $\sqrt{C_t}$)
- Transition function (next periods states)
- Value function = Value today (${V}^{optimal}_{t}({M}_{t})$)
- Continuation value (Value after today)
- Policy function (OPtimal choice)
