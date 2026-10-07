# Immune Model Scripts

## The Staged TIV model

The TIV model is a classic model of viral dynamics.  The staged TIV model extends the TIV framework to include staged progression in the I compartment.

$$\begin{align}
\frac{dT}{dt}    &= -\alpha T V \\
\frac{dV}{dt}    &= \sum_{i=m+1}^{m+n+1} p_i y_i - C V -\alpha T V\\
\frac{dy_1}{dt}  &= \alpha T V  - k_1 y_1 - d_1 y_1\\
\frac{dy_i}{dt}  &= k_{i-1} y_i - k_i y_i - d_i y_i \qquad i\in \{2,\dots,m+n+1\}
\end{align}$$

| state | symbol | notes |
|:------|:------:|:------|
| Target Cell | $T$ | Density of uninfected cells which are susceptible to infection |
| Virus | $V$ | Density of free virus, capable of infection |
| Infected Cell | $y_i$ | Density of infected cells in infection-stage $i$, $i\in\{1,\dots,m+n+1\}$ |

| parameter | symbol | notes |
|:------|:------:|:------|
| $\alpha$ |   |    |
| $p_i$ |   |    |
| $C$   |   |    |
| $k_i$ |   |    |
| $d_i$ |   |    |

## Formulation of the Stochastic Differential Equation Model


