## Installation

This package is not yet registered in the General registry. To install, type `]`
in the Julia command prompt, then execute

```julia
pkg> add https://github.com/ericphanson/SDPAFamily.jl
```

### Automatic binary installation

If you are on MacOS or Linux, this package will attempt to automatically
download the SDPA-GMP, SDPA-DD, and SDPA-QD binaries, built by
[SDPA\_GMP\_Builder.jl](https://github.com/ericphanson/SDPA_GMP_Builder). The
SDPA-GMP binary is patched from the official SDPA-GMP source to allow printing
more digits, in order to recover high-precision output.

SDPA-GMP does not compile on Windows. However, it can be used via the [Windows
Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about)(WSL).
If you have WSL installed, then SDPAFamily.jl will try to automatically detect
this and use an appropriate binary, called via WSL. This binary can be found at
the repo <https://github.com/ericphanson/SDPA_GMP_Builder>, and is built on WSL
from the source code at <https://github.com/ericphanson/sdpa-gmp>. Windows
support is experimental, however, and we could not get it to run on Travis. Any
help in this regard would be appreciated.

SDPA-{GMP, QD, DD} are each available under a GPLv2 license, which can be found
here:
<https://github.com/ericphanson/SDPA_GMP_Builder/blob/master/deps/COPYING>.

### Custom binary

If you would like to use a different binary, set the environmental variable
`JULIA_SDPA_GMP_PATH` (similarly for `JULIA_SDPA_QD_PATH` or
`JULIA_SDPA_DD_PATH`) to the folder containing the binary you would like to use,
and then build the package. This can be done in Julia by, e.g.,

```julia
ENV["JULIA_SDPA_GMP_PATH"] = "/path/to/folder/"
import Pkg
Pkg.build("SDPAFamily")
```

and that will configure this package to use that binary by default. If your
custom location is via WSL on Windows, then also set `ENV["JULIA_SDPA_GMP_WSL"]
= "TRUE"` (similarly for `ENV["JULIA_SDPA_QD_WSL"]` or
`ENV["JULIA_SDPA_DD_WSL"]`) so that SDPAFamily.jl knows to adjust the paths to
the right format. Note that the binary must be named `sdpa_gmp`, `sdpa_qd` and
`sdpa_dd`. 

It is recommended to patch SDPA-{GMP, QD, DD} (as was done in
<https://github.com/ericphanson/sdpa-gmp>) in order to allow printing more
digits. To do this for SDPA-GMP, and similarly for -QD and -DD,

* For source code downloaded from the official website (dated 20150320), modify
  the `P_FORMAT` string at line 23 in `sdpa_struct.h` so that the output has a
  precision no less than 200 bits (default) or precision specified by the
  parameter file. 
* For source code downloaded from its [GitHub
  repository](https://github.com/nakatamaho/sdpa-gmp), specify the print format
  string in `param.sdpa` as described in the [SDPA users
  manual](https://sourceforge.net/projects/sdpa/files/sdpa/sdpa.7.1.1.manual.20080618.pdf).

Other information about compiling SDPA solvers can be found
[here](http://sdpa.sourceforge.net/download.html). 
