## Non-Linear Microproblem

$$\frac{\partial c}{\partial t} = \nabla \cdot (\mathbf{D} \nabla c) - \frac{\partial f(c)}{\partial t}$$

Applying chain rule to $\frac{\partial f}{\partial t}$ yields

$$\frac{\partial c}{\partial t} = \nabla \cdot (\mathbf{D} \nabla c) - \frac{\partial f(c)}{\partial c}\frac{\partial c}{\partial t}.$$

Discretizing the time derivative by means of an implicit Euler scheme leads to

$$c^{n+1} - \Delta t \nabla \cdot (\mathbf{D} \nabla c^{n+1}) + \frac{\partial f(c^{n+1})}{\partial c^{n+1}} (c^{n+1}-c^{n}) - c^n = 0.$$

Integrating over domain $\Omega$ and multiplying test function $v$ results in

$$F:=\int_{\Omega} c^{n+1} \ v + \Delta t \ \nabla v^T \mathbf{D} \nabla c^{n+1} + {\color{red}\frac{\partial f(c^{n+1})}{\partial c^{n+1}}(c^{n+1}-c^n)\ v}-c^n\ v =0.$$

Using the following definitions
$$c= \sum_{j}^n \phi_i c_i$$
$$v = \phi_j \qquad \text{since } \forall \ v \in V_h $$
$$f:= \frac{Q_mK_nc}{1+K_n c}$$
$$\frac{\partial f(c)}{\partial c} = f'(c)=\frac{Q_m K_n}{(1+K_n c)^2} = Q_mK_n(1+K_nc)^{-2}$$
$$\frac{\partial^2 f(c)}{\partial c^2} =f''(c)= -2Q_mK_n^2(1+K_n\ c)^{-3}$$

we obtain the nonlinear problem with Jacobian


\begin{align}
J_{ij} := \frac{\partial F_i}{\partial c_i} = & \int_{\Omega} \phi_i\phi_j + \Delta t\ \nabla \phi_j ^T \mathbf{D} \nabla \phi_i + \\
 & {\color{red} f'(c^-) \phi_i \phi_j + f''(c^-)c^-\phi_i \phi_j - f''(c^-)c^n\phi_i\phi_j}
\end{align}

with $c^-$ being the most recent Newton approximation. 

$$\mathbf{J} \mathbf{\Delta c} = -\mathbf{F(c^-)}$$
$$\mathbf{c} = \mathbf{c}^- + \omega \mathbf{\Delta c}$$
$$\mathbf{c}^-\leftarrow \mathbf{c}$$
