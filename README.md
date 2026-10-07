# Immune Model Scripts

## The Staged TIV model

The TIV model is a classic model of viral dynamics.  The staged TIV model extends the TIV framework to include staged progression in the I compartment.

$$\begin{align}
\frac{dT}{dt}    &= -\alpha T V \\
\frac{dV}{dt}    &= \sum_{i=m+1}^{m+n+1} p_i y_i - C V -\alpha T V\\
\frac{dI_1}{dt}  &= \alpha T V  - k_1 I_1 - d_1 I_1\\
\frac{dI_i}{dt}  &= k_{i-1} I_i - k_i I_i - d_i I_i \qquad i\in \{2,\dots,m+n+1\}
\end{align}$$

The model breaks the infected compartment into $n+m+1$ compartments where $m$ is the number of eclipse stages and $n$ is the number of budding ramp-up stages.  Cells in the eclipse stages are infected, but not productively infected ($p_i =0$).  Cells in the ramp-up stages are productively infected, but not at maximum productivity ($p_i < p_{m+n+1}$.  Setting $m=1$, $n=0$ gives the classic TEIV model with a single eclipse stage and a single infectious stage.  Setting $m=n=0$ gives the classic TIV model.

| state | symbol | notes |
|:------|:------:|:------|
| Target Cell | $T$ | Density of uninfected cells which are susceptible to infection |
| Virus | $V$ | Density of free virus, capable of infection |
| Infected Cell | $I_i$ | Density of infected cells in infection-stage $i$, $i\in\{1,\dots,m+n+1\}$ |

| parameter | symbol | notes |
|:------|:------:|:------|
| infection rate       | $\alpha$ |   |
| budding rate         | $p_i$    |   |
| viral clearance rate | $C$      |   |
| progression rate     | $k_i$    |   |
| cell death rate      | $d_i$    |   |
| eclipse stages       | $m$      | $p_i=0$ for $i\in\{1,\dots,m\}$  |
| ramp-up stages       | $n$      | $p_{i}<p_{i+1}$ for $i\in\{m+1,\dots,m+n\}$ |

## Formulation of the Stochastic Differential Equation Model

The stochastic DE is formulated as

$$dx = f(x)dt + G(x)dW(t)$$

where $f$ is the vector field of the ODE model and $G$ is the square root of a covariance matrix.  In brief, $G$ has one column for each state transition event and the entries of $G$ consist of square roots of the transition rates.  Note the vector field consists of sums of the same transition rates.   For the staged TIV model with $m$ eclipse stages and $n$ ramp up stages there are $4 + 2(m+n)$ independent transition events.

| event           | transitions            |  rate        |   |
|:----------------|:----------------------:|:------------:|:--|
|infection        | $T-=1$; $V-=1$; $I_1+=1$   | $\alpha TV$   |
|budding          | $V+=1$                     | $\sum p_iI_i$    |
|free virus decay | $V-=1$                     | $C V$   |
|progression      | $I_i -= 1$; $I_{i+1} += 1$ | $k_i I_i$ | $i\in\{1,\dots,n+m\}$ |
|cell death       | $I_i -=1$                  | $d_i I_i$ | $i\in\{1,\dots,n+m+1\}$   |

Each column of $G$ corresponds to a state variable, and each row of $G$ corresponds to a transition event.  If we let $\zeta$ denote the matrix whose columns denote the transitions by $\pm 1$, then $G = \zeta \text{diag}(\sqrt{f}}$ so that $GG^T$ is the (square) covariance matrix.
