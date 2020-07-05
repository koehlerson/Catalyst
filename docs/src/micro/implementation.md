## Structs
Suitable datastructures are developed that hold as their fields all information about the microstructure. 
In case of the simplified, linerized sink term `CatalystStateODE` can be used to store various numbers associated with the catalyst state. 
When a PDE description of the microscale is used `CatalystStatePDE` should be used. This struct assembles by a given micro diffusion the corresponding linear operators of the micro scale and optionally assembles in each `catalyst_update!` call the nonlinear parts if `:nonlinear` is passed as the `microcomp_type`. 

```@docs
CatalystStateODE
CatalystStatePDE
```

## Assembly

Most of the assembly is handled while constructing the corresponding CatalystState struct. However, in case of a nonlinear computation the nonlinear parts, i.e. the Jacobian of the residual has to be assembled in each newton iteration, due to its nonlinear nature.

```@docs
Catalyst.assemble_nonlinear_micro_global!
```
!!! note "Care the interface"
    `Catalyst.assemble_nonlinear_micro_global!` is the only assembly routine which mutates the input arguments Jacobi $\mathbf{K}$ and residual $\mathbf{f}$
```@docs
Catalyst.assemble_nonlinear_micro_element!
```
## Micro Solver 
From a user perspective only `catalyst_update` is invoked which determines based on the input which micro computation need to be done. In case of an array of arrays holding `CatalystStateODE`s the update is done based on a suitable ODE. 

If a array of array of `CatalystStatePDE` is passed, either a linear or nonlinear micro PDE is solved.

```@docs
Catalyst.catalyst_update!
```
### Linear Micro Solver
```@docs
Catalyst.microcomputation_linear!
```
### Nonlinear Micro Solver
```@docs
Catalyst.microcomputation_nonlinear!
```

## Micro Source/Sink

```@docs
Catalyst.langmuir_isotherm′
Catalyst.langmuir_isotherm″
```
