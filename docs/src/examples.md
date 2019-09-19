# Examples

```@setup 1
# This setup block is not shown in the final output
# Install the right branch of Convex
using Pkg
Pkg.add(PackageSpec(name="Convex", url="https://github.com/ericphanson/Convex.jl", rev="MathOptInterface"));
```

Here is a simple optimization problem formulated with Convex.jl:

```@example 1
using SDPAFamily, LinearAlgebra
using Convex # ] add https://github.com/ericphanson/Convex.jl#MathOptInterface
y = Semidefinite(3)
p = maximize(lambdamin(y), tr(y) <= 5; numeric_type = BigFloat)
solve!(p, SDPAFamily.Optimizer(presolve=true))
@show p.optval
```

## Optimal guessing probability for a pair of quantum states

In physics, a *state* represents a possible configuration of a physical system.
In quantum mechanical systems with finitely many degrees of freedom, states are
represented by *density matrices*, which are $d\times d$ matrices with complex
entries that are positive semi-definite and have trace equal to one. States can
be *measured*; mathematically, a measurement with $n$ possible outcomes is
represented by a set of measurement operators $\{E_j\}_{j=1}^n$, where each
$E_j$ is a $d\times d$ matrix. For example, imagine an experiment in which a
charged particle is released in a magnetic field such that it will hit either a
detector on the left or a detector on the right. This corresponds to a
measurement of the particle with two outcomes, and hence two measurement
operators $\{E_1, E_2\}$, which to the left and right detector.

In order for $\{E_j\}_{j=1}^n$ to be a valid set of measurement operators, each
$E_j$ must be positive semi-definite, and the family $\{E_j\}_{j=1}^n$ must have
the property that $\sum_{j=1}^n E_j = I_d$, the $d\times d$ identity matrix. If
the state of the system is represented by $\rho$, and a measurement with
measurement operators $\{E_j\}_{j=1}^n$ is performed, then outcome $j$ is
obtained with probability $\operatorname{tr}[\rho E_j]$.

Consider the case where $d=2$ (i.e. the states are *qubits*), and the state of
the system is either represented by $\rho_1 = \begin{pmatrix} 1 & 0 \\ 0 & 0
\end{pmatrix}$ or by $\rho_2 = \frac{1}{2}\begin{pmatrix} 1 & -i \\ i & 1
\end{pmatrix}$, but we don't know which; let's say there is a referee who
flipped a fair coin, and then prepared the system in either $\rho_1$ or
$\rho_2$. We will perform a measurement of the system, and then use the outcome
to make a guess about the state of the system.  What is the measurement that
gives the highest probability of correctly determining which state the system is
in, and what's the optimal probability?

We will perform a measurement with measurement operators $E_1$ and $E_2$. If we
get outcome $1$, we will guess the system is in state $\rho_1$ and and if we
obtain outcome 2, we guess the system is in state $\rho_2$. Then the probability
of guessing correctly is

```math
p_\text{guess}(E_1, E_2) = \frac{1}{2}\operatorname{tr}(\rho_1  E_1) + \frac{1}{2}\operatorname{tr}(\rho_2  E_2)
```

since there is a 50% chance of the system being in state $\rho_1$, in which case
we guess correctly when we get outcome 1 (which occurs with probability
$\operatorname{tr}(\rho_1 E_1)$), and a 50% chance of the system being in state
$\rho_2$, in which case we guess correctly when we get outcome $2$.

Our goal now is to choose the optimal measurement operators to have the the best
chance of guessing correctly. That is, we aim to maximize the above expression
over all choices of $E_1$ and $E_2$ such that $\{E_1, E_2\}$ is a valid set of
measurement operators. This is a semidefinite program, which can be solved e.g.
with SDPAFamily.jl In this simple example with only two states to discriminate
between, the problem can be solved analytically, and the solution is related to
the trace distance between the two states. This problem specifically is Example
3.2.1 of the [edX Quantum Cryptography notes by Thomas
Vidick](http://users.cms.caltech.edu/~vidick/teaching/120_qcrypto/LN_Week3.pdf).
It can be seen that the optimal guessing probability is

```math
p_\text{guess} = \frac{1}{2} + \frac{1}{2 \sqrt{2}}
```

Let us see to what accuracy we can recover that result using the SDPA solvers.

```@example 1
using SDPAFamily, Printf
using Convex # ] add https://github.com/ericphanson/Convex.jl#MathOptInterface

ρ₁ = Complex{BigFloat}[1 0; 0 0]
ρ₂ = (1//2)*Complex{BigFloat}[1 -im; im 1]

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

Here, we have solved the problem four times, once with each variant of the SDPA
family of optimizers supported by this package. We can see that SDPA-GMP has
solved the problem to an accuracy of $\sim 10^{-31}$, far exceeding machine
precision.

As usual with semidefinite programs, we can recover a set of optimal
measurements:

```@example 1
evaluate(E₁)
```

```@example 1
evaluate(E₂)
```

Note that this is an example where the presolve routine is essential to getting
good results:

```@example 1
for variant in (:sdpa, :sdpa_dd, :sdpa_qd, :sdpa_gmp)
    solve!(problem, SDPAFamily.Optimizer(silent = true, presolve = false, variant = variant))
    error = abs(problem.optval - p_guess)
    print("$variant solved the problem with an absolute error of ")
    @printf("%.2e.\n", error)
end
```

We can see that without the presolve routine, we have only recovered the true
solution up to errors of size $\sim 10^{-1}$ for `:sdpa` variant. All other
variants have failed to produce a result due to redundant constraints and
returned with default value 0.
