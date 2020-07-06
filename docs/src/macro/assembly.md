This page covers assembly routines that are used for the macroscopic problem. 
The interface of all assembly routines is the same in that sense that, depending on the input arguments (multiple dispatch) doassemble returns a discrete FE Operator which can be used to set up the linear system of equation or the differential algebraic equation in time, respectively.

More about the usage with the discrete operators can be found in the solver section.

```@docs
doassemble
Catalyst.volume
```
