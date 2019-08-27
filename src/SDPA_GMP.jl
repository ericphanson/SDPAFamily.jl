module SDPA_GMP
using LinearAlgebra # for diag()
using DelimitedFiles # for writing SDPA input file
using MathOptInterface
MOI = MathOptInterface
MOIB = MOI.Bridges
# const sdpa_gmp_path = "sdpa_gmp" # make sure sdpa_gmp binary is in the system PATH

export sdpa_gmp_binary_solve!, read_results!

include("MOI_wrapper.jl")
include("file_io.jl")
include("presolve.jl")
include("binary_call.jl")

end # module
