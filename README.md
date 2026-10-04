# Lorenz Attractor: A Mathematical Journey

This project makes an interactive 3D picture of the Lorenz system. I started
with a question that looks simple: how can a rule that always gives the same
answer create motion that is difficult to predict far ahead?

The graph is the system's path through **phase space**. Each point is the
system's state at one moment, and the three coordinates change according to
the same differential equations throughout the run. The butterfly shape is
not a butterfly moving through ordinary space.

## Why the butterfly is interesting

The path circles around one wing, switches to the other, and keeps doing this
without repeating the same loop. Starting points that are almost identical
can eventually lead to noticeably different paths. The system is still
deterministic: if I knew the exact starting point and had exact arithmetic, the
rule would specify the trajectory. In practice, tiny differences grow, making
long-term prediction difficult.

This is a simplified mathematical model associated with convection, not a
complete weather simulation. The colorful curve shows time: hovering over it
reveals the time and the corresponding \((x,y,z)\) coordinates. The white
markers show equilibria, with Jacobian eigenvalues available on hover.

## The equations

The state is \(\gamma(t)=(x(t),y(t),z(t))\). Its motion is set by the Lorenz
vector field:

$$
\begin{aligned}
\dot{x} &= \sigma(y-x),\\
\dot{y} &= x(\rho-z)-y,\\
\dot{z} &= xy-\beta z.
\end{aligned}
$$

For this picture, the parameters are \(\sigma=10\), \(\rho=28\), and
\(\beta=8/3\). The path starts at \((1,1,1)\). A dot means a derivative with
respect to time, so \(\dot{x}\), for example, is the rate at which \(x\) is
changing at that instant.

## My analysis walkthrough

### 1. Read the equations as a direction field

At every point \((x,y,z)\), the equations give a velocity vector
\(F(x,y,z)\). The trajectory follows that vector:

$$
\gamma'(t)=F(\gamma(t)).
$$

Writing each equation with differentials makes the small changes explicit:

$$
\begin{aligned}
dx &= \sigma(y-x)\,dt,\\
dy &= \bigl(x(\rho-z)-y\bigr)\,dt,\\
dz &= (xy-\beta z)\,dt.
\end{aligned}
$$

My first expectation is that one exact starting state gives one path. The
interesting question is what happens to a *nearby* starting state.

### 2. Try a curve integral

I can turn the vector field into a one-form and integrate it along the path:

$$
\omega=F_x\,dx+F_y\,dy+F_z\,dz.
$$

Because the path's tangent is the vector field itself,
\(\gamma'(t)=F(\gamma(t))\), this particular integral becomes the accumulated
squared speed in these coordinates:

$$
\begin{aligned}
\int_\gamma \omega
&=\int_0^{60}F(\gamma(t))\cdot\gamma'(t)\,dt\\
&=\int_0^{60}\|F(\gamma(t))\|^2\,dt\\
&\approx 6.2589\times 10^5.
\end{aligned}
$$

The number is calculated from the sampled numerical trajectory. It is a useful
worked example of a curve integral, but it is **not** physical work or energy:
the coordinates describe the model's state, and changing their scale changes
the value.

### 3. Ask whether a little cloud expands or shrinks

I imagine a tiny 3D cloud of starting states around one point. The divergence
summarizes the flow's instantaneous net effect on the cloud's volume. For this
system it is constant:

$$
\nabla\cdot F=-(\sigma+1+\beta)=-\frac{41}{3}.
$$

If \(V(t)\) is a sufficiently small volume carried by the flow, then

$$
\frac{dV}{dt}=(\nabla\cdot F)V,
\qquad
\frac{V(t+1)}{V(t)}=e^{-41/3}\approx1.16\times10^{-6}.
$$

So nearby states are squeezed together **in total volume**. That does not mean
every direction shrinks or every pair of trajectories gets closer. Some
directions can stretch while other directions contract.

### 4. Use the Jacobian to look locally

The Jacobian records how the vector field changes when I make a small change
in the state. A tiny displacement \(\delta\) evolves approximately as

$$
\delta'(t)=DF(\gamma(t))\,\delta(t),
\qquad
DF(x,y,z)=
\begin{pmatrix}
-\sigma & \sigma & 0\\
\rho-z & -1 & -x\\
y & x & -\beta
\end{pmatrix}.
$$

This is a local linear approximation. It helps explain stretching and
shrinking near the current point, but the matrix changes as the trajectory
moves. The combination of local stretching, contraction, and folding helps
explain the complicated bounded shape. Negative divergence tells me about
volume contraction; on its own, it does not prove chaos.

## Quick run

From the project folder, install the Python packages and run the script:

```powershell
python -m pip install numpy plotly
python lorenz_attractor.py
```

The script creates `lorenz_attractor.html` and attempts to open it in your
default browser. The Plotly graph is interactive: rotate, zoom, pan, inspect
coordinates and time, and hover over the equilibria. The formulas are rendered
with KaTeX from a CDN, so an internet connection is needed to display them.

If `python` points to a different environment than the one where the packages
were installed, use that environment's Python executable for both commands.

## Project files

```text
LORENZ/
├── lorenz_attractor.py   # Lorenz equations, RK4 integration, HTML generation
├── lorenz_template.html  # Page layout, explanations, and KaTeX formulas
├── lorenz_attractor.html # Generated interactive graph and explanations
└── README.md
```

## What the code currently does

The script uses the classic parameter values \((\sigma,\rho,\beta)=(10,28,8/3)\),
starts at \((1,1,1)\), and takes 12,000 fourth-order Runge-Kutta steps with
\(dt=0.005\). This covers time from \(t=0\) to \(t=60\). It calculates the
equilibria and their Jacobian eigenvalues for the hover details, and computes
the curve-integral example from the same sampled path.

## Interpretation note

This is a numerical illustration of one trajectory, not a proof that every
starting point behaves the same way. The familiar butterfly depends on the
parameter values, initial condition, and numerical approximation. The
divergence result is exact for these parameters; the displayed curve integral
is a numerical approximation for this particular run.
