# Examples

## Optimal guessing probability for a pair of quantum states

```@setup 1
# This setup block is not shown in the final output
# Install the right branch of Convex
using Pkg
Pkg.add(PackageSpec(name="Convex", url="https://github.com/ericphanson/Convex.jl", rev="MathOptInterface"));
```

In physics, a *state* represents a possible configuration of a physical system. In quantum mechanical systems with finitely many degrees of freedom, states are represented by *density matrices*, which are $d\times d$ matrices with complex entries that are positive semi-definite and have trace equal to one. States can be *measured*; mathematically, a measurement with $n$ possible outcomes is represented by a set of measurement operators $\{E_j\}_{j=1}^n$. Each $E_j$ is a $d\times d$ matrix which is positive semi-definite, and the family $\{E_j\}_{j=1}^n$ has the property that $\sum_{j=1}^n E_j = I_d$, the $d\times d$ identity matrix. If the state of the system is represented by $\rho$, and a measurement with measurement operators $\{E_j\}_{j=1}^n$ is performed, then outcome $j$ is obtained with probability $\operatorname{tr}[\rho E_j]$.

Consider the case where $d=2$ (i.e. the states are *qubits*), and the state of the system is either represented by $\rho_1 = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$ or by $\rho_2 = \frac{1}{2}\begin{pmatrix} 1 & -i \\ i & 1 \end{pmatrix}$, but we don't know which. We have a prior assumption that there is a 50% chance the state of the system is given by $\rho_1$ and a 50% chance the state of the system is given by $\rho_2$. What is the measurement that gives the highest probability of correctly determining which state the system is in?

Assume the measurement operators are $E_1$ and $E_2$, and if we obtain outcome 1, we guess the system is in state $\rho_1$ and and if we obtain outcome 2, we guess the system is in state $\rho_2$, then the probability of guessing correctly is

```math
\frac{1}{2}\operatorname{tr}(\rho_1  E_1) + \frac{1}{2}\operatorname{tr}(\rho_2  E_2)
```

since there is a 50% chance of the system being in state $\rho_1$, in which case we guess correctly when we get outcome 1 (which occurs with probability $\operatorname{tr}(\rho_1 E_1)$), and a 50% chance of the system being in state $\rho_2$, in which case we guess correctly when we get outcome $2$.

Our goal now is to choose the optimal measurement operators to have the the best chance of guessing correctly. That is, we aim to maximize the above expression over all choices of $E_1$ and $E_2$ such that $\{E_1, E_2\}$ is a valid set of measurement operators. This is a semidefinite program, which can be solved e.g. with SDPAFamily.jl In this simple example, the problem can be solved analytically; in fact, this problem is Example 3.2.1 of the [edX Quantum Cryptography notes by Thomas Vidick](http://users.cms.caltech.edu/~vidick/teaching/120_qcrypto/LN_Week3.pdf), from which it can be seen that the correct answer is

```math
p_\text{guess} = \frac{1}{2} + \frac{1}{2 \sqrt{2}}
```

Let us see to what accuracy we can recover that result using the SDPA solvers.

```@example 1
using SDPAFamily, Printf
using Convex # ] add https://github.com/ericphanson/Convex.jl#MathOptInterface

ρ₁ = [big"1.0" big"0.0"; big"0.0" big"0.0"]
ρ₂ = [big"0.5" big"-0.5"*im; big"0.5"*im big"0.5"]

E₁ = ComplexVariable(2, 2)
E₂ = ComplexVariable(2, 2)

problem = maximize(real((1//2)*tr(ρ₁*E₁) + (1//2)*tr(ρ₂*E₂)),
            [E₁ ⪰ 0, E₂ ⪰ 0, E₁ + E₂ == Diagonal(ones(2))]; numeric_type = BigFloat)
p_guess = 1//2 + 1/(2*sqrt(big(2)))
for variant in (:sdpa, :sdpa_dd, :sdpa_qd, :sdpa_gmp)
    solve!(problem, SDPAFamily.Optimizer(silent = true, presolve = true, variant = variant))
    error = abs(problem.optval - p_guess)
    print("$variant solved the problem with an absolute error of ")
    @printf("%.2e.\n", error)
end
```

As usual with semidefinite programs, we can recover a set of optimal measurements:

```@example 1
evaluate(E₁)
```

```@example 1
evaluate(E₂)
```

Note that this is an example where the presolve routine is essential to getting good results:

```@example 1
for variant in (:sdpa, :sdpa_dd, :sdpa_qd, :sdpa_gmp)
    solve!(problem, SDPAFamily.Optimizer(silent = true, presolve = false, variant = variant))
    error = abs(problem.optval - p_guess)
    print("$variant solved the problem with an absolute error of ")
    @printf("%.2e.\n", error)
end
```
