# Immune Model Scripts

## The Staged TRIV model

The TIV model is a classic model of viral dynamics.  The staged TIV model extends the TIV framework to include staged progression in the I compartment.
The TRIV model includes a compartment of *refractory* cells.   There are many variations on the base TRIV model:  [@ciupe2007Role] assumes a flow from I to R in proportion to a population of immune effector cells;  whereas [@ke2022Daily] assumes a flow from T to R in proportion to $IT$.   Here, for simplicity, we assume cells move from $T$ to $R$ in proportion to the viral load.  This glosses over the details of cytokine signalling leading to resistance.

$$\begin{align}
\frac{dT}{dt}    &= -\alpha T V - \phi T V \\
\frac{dR}{dt}    &= \phi T V - \rho R\\
\frac{dV}{dt}    &= \sum_{i=m+1}^{m+n} p_i y_i - C V -\alpha T V\\
\frac{dy_1}{dt}  &= \alpha T V  - k_1 y_1 - d_1 y_1\\
\frac{dy_i}{dt}  &= k_{i-1} y_i - k_i y_i - d_i y_i \qquad i\in \{2,\dots,m+n\}
\end{align}$$


