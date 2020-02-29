# Examples

```@setup 1
# This setup block is not shown in the final output
setprecision(256)
```

Here is a simple optimization problem formulated with Convex.jl:

```@repl 1
using SDPAFamily, LinearAlgebra
using Convex
y = Semidefinite(3)
p = maximize(lambdamin(y), tr(y) <= 5; numeric_type = BigFloat)
solve!(p, () -> SDPAFamily.Optimizer(presolve=true))
p.optval
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

```@repl 1
using Convex, SDPAFamily, Printf

ρ₁ = Complex{BigFloat}[1 0; 0 0]
ρ₂ = (1//2)*Complex{BigFloat}[1 -im; im 1]
E₁ = ComplexVariable(2, 2)
E₂ = ComplexVariable(2, 2)
problem = maximize( real((1//2)*tr(ρ₁*E₁) + (1//2)*tr(ρ₂*E₂)),
                    [E₁ ⪰ 0, E₂ ⪰ 0, E₁ + E₂ == Diagonal(ones(2))];
                    numeric_type = BigFloat );
p_guess = 1//2 + 1/(2*sqrt(big(2)))
for variant in (:sdpa, :sdpa_dd, :sdpa_qd, :sdpa_gmp)
    solve!(problem, () -> SDPAFamily.Optimizer(silent = true, presolve = true, variant = variant))
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

```@repl 1
evaluate(E₁)
evaluate(E₂)
```

Note that this is an example where the presolve routine is essential to getting
good results:

```@repl 1
for variant in (:sdpa, :sdpa_dd, :sdpa_qd, :sdpa_gmp)
    solve!(problem, () -> SDPAFamily.Optimizer(silent = true, presolve = false, variant = variant))
    error = abs(problem.optval - p_guess)
    print("$variant solved the problem with an absolute error of ")
    @printf("%.2e.\n", error)
end
```

We can see that without the presolve routine, we have only recovered the true
solution up to errors of size $\sim 10^{-1}$ for `:sdpa` variant. All other
variants have failed to produce a result due to redundant constraints and
returned with default value 0.

This problem is revisited at very high precision in [Changing parameters & solving at very high precision](@ref).

## Polynomial optimization

The following example is adapted from an example in the
[SumOfSquares.jl](https://github.com/JuliaOpt/SumOfSquares.jl) documentation to
use SDPAFamily.jl. Even though the problem is only specified with `Float64`'s, since
the entries are specified as integers, they can be sent to SDPA-GMP without a loss of
precision.

```@repl sumofsquares
using SumOfSquares
using DynamicPolynomials
using SDPAFamily

@polyvar x1 x2 # Create symbolic variables (not JuMP decision variables)

# Create a Sum of Squares JuMP model with the SDPAFamily solver
model = SOSModel(with_optimizer(SDPAFamily.Optimizer{Float64}, # JuMP only supports Float64
                    variant = :sdpa_gmp, # use the arbitrary precision variant
                    params = (  epsilonStar = 1e-30, # constraint tolerance
                                epsilonDash = 1e-30, # normalized duality gap tolerance
                                precision = 200 # arithmetric precision used in sdpa_gmp
                )))

@variable(model, γ) # Create a JuMP decision variable for the lower bound

# f(x) is the Goldstein-Price function
f1 = x1 + x2 + 1
f2 = 19 - 14 * x1 + 3 * x1^2 - 14 * x2 + 6 * x1 * x2 + 3 * x2^2
f3 = 2 * x1 - 3 * x2
f4 = 18 - 32 * x1 + 12 * x1^2 + 48 * x2 - 36 * x1 * x2 + 27 * x2^2

f = (1 + f1^2 * f2) * (30 + f3^2 * f4)

@constraint(model, f >= γ) # Constrains f(x) - γ to be sum of squares

@objective(model, Max, γ)

optimize!(model)

println(objective_value(model))
```

Let's check the input file that is used to specify the problem for the SDPA-GMP
binary. The following command uses some implementation details of how JuMP
stores the underlying optimizer and so may not work in later JuMP versions.
However, the following path is always printed when setting `verbose =
SDPAFamily.VERBOSE`.

```@repl sumofsquares
path = joinpath(MOI.get(model, SDPAFamily.TemporaryDirectory()),  "input.dat-s")
readlines(path)[1:10] .|> println;
```

The full file is longer, but what gets passed to the optimizer for this problem
are floating point numbers that can be faithfully read by SDPA-GMP at the
200-bits of precision it uses internally. Thus, in this case, that JuMP
restricts the Julia model to store the numbers at machine precision does not
affect the precision of data that SDPA-GMP receives. It does, however, affect
the precision of data that JuMP can recover from the output file. In this case,
JuMP receives the correct answer to full machine precision ($\sim 10^{-16}$),
but the true answer printed by SDPA-GMP (which can be seen in the file
`output.dat-s`) is in fact correct to $\sim 10^{-30}$ in this case.

For this kind of problem which uses JuMP, the precision advantage of SDPA-GMP
over other problems is that SDPA-GMP should be able to solve the problem to the
full $\sim 10^{-16}$ precision representable by 64-bit floating point numbers,
while solvers which solve the problem in machine precision can only recover the
result to $\sim 10^{-8}$.

## SDP relaxation in polynomial optimization problem

We consider a polynomial optimization problem (POP) where we wish to find $p^* = \inf \lbrace x | x \geq 0,\ x^2 \geq 1 \rbrace$, which has an SDP relaxation of order $r$ as

```@raw html
<img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\begin{align*}&space;\text{minimize}&\&space;y_1&space;\\&space;\text{such&space;that}&&space;\begin{bmatrix}&space;1&space;&y_1&space;&\dots&space;&y_r&space;\\&space;y_1&space;&&space;y_2&space;&\dots&space;&y_{r&plus;1}&space;\\&space;\vdots&space;&&space;\vdots&space;&&space;\ddots&space;&\vdots&space;\\&space;y_r&space;&&space;y_{r&plus;1}&space;&\dots&space;&y_{2r}&space;\end{bmatrix}&space;&&space;\succeq&space;O\\&space;&&space;\begin{bmatrix}&space;y_1&space;&y_2&space;&\dots&space;&y_r&space;\\&space;y_2&space;&&space;y_3&space;&\dots&space;&y_{r&plus;1}&space;\\&space;\vdots&space;&&space;\vdots&space;&&space;\ddots&space;&\vdots&space;\\&space;y_r&space;&&space;y_{r&plus;1}&space;&\dots&space;&y_{2r-1}&space;\end{bmatrix}&space;&&space;\succeq&space;O\\&space;&\begin{bmatrix}&space;y_2&space;-&space;1&space;&y_3&space;-&space;y_1&space;&\dots&space;&y_{r&plus;1}&space;-&space;y_{r-1}&space;\\&space;y_3-y_1&space;&&space;y_4-y_2&space;&\dots&space;&y_{r&plus;2}&space;-&space;y_r&space;\\&space;\vdots&space;&&space;\vdots&space;&&space;\ddots&space;&\vdots&space;\\&space;y_{r&plus;1}&space;-&space;y_{r-1}&space;&&space;y_{r&plus;2}&space;-y_{r}&space;&\dots&space;&y_{2r}&space;-&space;y_{2r-2}&space;\end{bmatrix}&space;&&space;\succeq&space;O\\&space;&(y_1,\dots,y_{2r})&space;\in&space;\mathbb{R}^{2r}&space;\end{align*}" title="\large \begin{align*} \text{minimize}&\ y_1 \\ \text{such that}& \begin{bmatrix} 1 &y_1 &\dots &y_r \\ y_1 & y_2 &\dots &y_{r+1} \\ \vdots & \vdots & \ddots &\vdots \\ y_r & y_{r+1} &\dots &y_{2r} \end{bmatrix} & \succeq O\\ & \begin{bmatrix} y_1 &y_2 &\dots &y_r \\ y_2 & y_3 &\dots &y_{r+1} \\ \vdots & \vdots & \ddots &\vdots \\ y_r & y_{r+1} &\dots &y_{2r-1} \end{bmatrix} & \succeq O\\ &\begin{bmatrix} y_2 - 1 &y_3 - y_1 &\dots &y_{r+1} - y_{r-1} \\ y_3-y_1 & y_4-y_2 &\dots &y_{r+2} - y_r \\ \vdots & \vdots & \ddots &\vdots \\ y_{r+1} - y_{r-1} & y_{r+2} -y_{r} &\dots &y_{2r} - y_{2r-2} \end{bmatrix} & \succeq O\\ &(y_1,\dots,y_{2r}) \in \mathbb{R}^{2r} \end{align*}" align="middle"/>
```

This is an example when most solvers fail due to numerical errors using Float64. It can be shown that for all $r \geq 1$, we have $y_1^* = 0$ as the optimal value for the SDP relaxation. However, most solvers will report  $y_1^* = 1$, which is in fact the optimal value for $p^*$ in the original problem. The details are discussed in [1]. Such problem can be overcome by using appropriate parameters for `SDPA-GMP`. We now demonstrate this using `Convex.jl`.

```@repl 1
using SDPAFamily, SCS, Convex

function relaxed_pop(r::Int, T)
    v = Variable(2*r)
    M1 = v[1:1+r]'
    for i in 2:r
        M1 = vcat(M1, v[i:i+r]')
    end
    t = [1 v[1:r]']
    M1 = vcat(t, M1)
    c1 = M1 in :SDP
    M2 = M1[2:end, 1:end-1]
    c2 = M2 in :SDP
    M3 = M1[2:end, 2:end] - M1[1:end-1, 1:end-1]
    c3 = M3 in :SDP
    return Problem{T}(:minimize, v[1], [c1, c2, c3])
end

p1 = relaxed_pop(5, Float64);
solve!(p1, () -> SCS.Optimizer(max_iters = 10000, verbose = 0));
p2 = relaxed_pop(5, BigFloat);
solve!(p2, () -> SDPAFamily.Optimizer(presolve = true, verbose = SDPAFamily.SILENT,
            params = ( epsilonStar = 1e-90,
                       epsilonDash = 1e-90,
                       precision = 5000,
                       betaStar = 0.5,
                       betaBar = 0.5,
                       gammaStar = 0.5,
                       lambdaStar = 1e5,
                       omegaStar = 2.0,
                       maxIteration = 10000
                       )));
                       
p1.status
p1.optval
p2.status
p2.optval
```

[1] H. Waki, M. Nakata, and M. Muramatsu, ‘Strange behaviors of interior-point methods for solving semidefinite programming problems in polynomial optimization’, *Comput Optim Appl*, vol. 53, no. 3, pp. 823–844, Dec. 2012.
