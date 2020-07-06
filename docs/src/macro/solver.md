## Theory 

After assembling the discrete FE operators we are left with large system of ODEs in time, i.e.

$$\mathbf{M} \mathbf{c}'(t) + \mathbf{A} \mathbf{c}(t) = \mathbf{r}(\mathbf{c})$$

where $\mathbf{M}, \mathbf{A}$ and $\mathbf{r}$ represents the assembled mass matrix, diffusion & advection matrices and reaction operator, respectively. 
In this represenation it should become clear why we are left with a system of ODEs.
Furthermore note that the reaction operator is dependent of the current solution of $c$ and thus introduces a non-linearity.

This non-linearity is right now bypassed by an operator splitting and in this particular case by an sequential solving. So, we first solve the nonlinear micro and then use it and solve the linear problem.

## Implementation

```@docs
Catalyst.solve
```
