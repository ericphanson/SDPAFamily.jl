# Developer reference

## High level picture

The high precision SDPA solvers do not provide a library interface, only binary access. So SDPAFamily.jl takes the problem inputs from a MathOptInterface Optimizer object and writes them to a file in a temporary folder (in the SDPA format) and calls the binary. The binary reads the input file, solves the problem, and writes an output file. SDPAFamily.jl then reads the output file and populates the MathOptInterface Optimizer object.

## Docstrings

```@autodocs
Modules = [SDPAFamily]
Order   = [:type, :constant, :function]
```
