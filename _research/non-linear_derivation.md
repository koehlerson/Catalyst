## Non-Linear Microproblem

$$\frac{\partial c}{\partial t} = \nabla \cdot (\mathbf{D} \nabla c) - \frac{\partial f(c)}{\partial t}$$

Applying chain rule to $\frac{\partial f}{\partial t}$ yields

$$\frac{\partial c}{\partial t} = \nabla \cdot (\mathbf{D} \nabla c) - \frac{\partial f(c)}{\partial c}\frac{\partial c}{\partial t}.$$

Discretizing the time derivative by means of an implicit Euler scheme leads to

$$c^{n+1} - \Delta t \nabla \cdot (\mathbf{D} \nabla c^{n+1}) + \frac{\partial f(c^{n+1})}{\partial c^{n+1}} (c^{n+1}-c^{n}) - c^n = 0.$$

Integrating over domain $\Omega$ and multiplying test function $v$ results in

$$F:=\int_{\Omega} c^{n+1} \ v + \Delta t \ \nabla v^T \mathbf{D} \nabla c^{n+1} + \red{\frac{\partial f(c^{n+1})}{\partial c^{n+1}}(c^{n+1}-c^n)\ v}-c^n\ v =0.$$

Using the following definitions
$$c= \sum_{j}^n \phi_i c_i$$
$$v = \phi_j \qquad \text{since } \forall \ v \in V_h $$
$$f:= \frac{Q_mK_nc}{1+K_n c}$$
$$\frac{\partial f(c)}{\partial c} = \frac{Q_m K_n}{(1+K_n c)^2} = Q_mK_n(1+K_nc)^{-2}$$

we obtain the nonlinear problem with Jacobian


$$J_{ij} := \frac{\partial F_i}{\partial c_i^{n+1}} = \phi_i\phi_j + \Delta t\ \nabla \phi_j ^T \mathbf{D} \nabla \phi_i \red{ + Q_m K_n (1+K_nc^{n+1})^{-2}\phi_i \phi_j -2Q_mK_n^2(1+K_n\ c^{n+1})^{-3}c^{n+1}\phi_i\phi_j}$$
