import numpy as np
import plotly.graph_objects as go
from pathlib import Path


def lorenz_vector_field(state, sigma=10.0, rho=28.0, beta=8.0 / 3.0):
    """Return the Lorenz system's velocity at a point (x, y, z)."""
    x, y, z = state
    return np.array(
        [sigma * (y - x), x * (rho - z) - y, x * y - beta * z],
        dtype=float,
    )


def integrate_lorenz(initial_state, dt=0.005, steps=12000):
    """Approximate one orbit with the fourth-order Runge-Kutta method."""
    trajectory = np.empty((steps + 1, 3), dtype=float)
    trajectory[0] = initial_state

    for index in range(steps):
        state = trajectory[index]
        k1 = lorenz_vector_field(state)
        k2 = lorenz_vector_field(state + 0.5 * dt * k1)
        k3 = lorenz_vector_field(state + 0.5 * dt * k2)
        k4 = lorenz_vector_field(state + dt * k3)
        weighted_slopes = k1 + 2 * k2 + 2 * k3 + k4
        trajectory[index + 1] = state + (dt / 6.0) * weighted_slopes

    return trajectory


def plot_attractor(trajectory, dt, initial_state):
    """Save and open an interactive 3D plot coloured by time."""
    time = np.arange(len(trajectory)) * dt
    sigma, rho, beta = 10.0, 28.0, 8.0 / 3.0
    velocity = np.column_stack(
        (
            sigma * (trajectory[:, 1] - trajectory[:, 0]),
            trajectory[:, 0] * (rho - trajectory[:, 2]) - trajectory[:, 1],
            trajectory[:, 0] * trajectory[:, 1] - beta * trajectory[:, 2],
        )
    )
    speed_squared = np.einsum("ij,ij->i", velocity, velocity)
    curve_integral = np.sum(
        0.5 * (speed_squared[:-1] + speed_squared[1:]) * dt
    )
    integral_mantissa, integral_exponent = f"{curve_integral:.4e}".split("e")
    curve_integral_latex = (
        f"{integral_mantissa}\\times10^{{{int(integral_exponent)}}}"
    )

    figure = go.Figure()
    figure.add_trace(
        go.Scatter3d(
            x=trajectory[:, 0],
            y=trajectory[:, 1],
            z=trajectory[:, 2],
            mode="lines",
            name="Trajectory",
            line={
                "color": time,
                "colorscale": "Turbo",
                "width": 4,
                "colorbar": {"title": "Time"},
            },
            customdata=time,
            hovertemplate=(
                "t = %{customdata:.3f}<br>"
                "x = %{x:.4f}<br>y = %{y:.4f}<br>z = %{z:.4f}"
                "<extra>Trajectory</extra>"
            ),
        )
    )

    # I used the fixed points as reference spots to compare the orbit with the
    # local linear behaviour predicted by the Jacobian at each equilibrium.
    equilibria = [
        ("Origin", 0.0, 0.0, 0.0),
    ]
    wing_coordinate = np.sqrt(beta * (rho - 1.0))
    equilibria.extend(
        [
            ("C+", wing_coordinate, wing_coordinate, rho - 1.0),
            ("C-", -wing_coordinate, -wing_coordinate, rho - 1.0),
        ]
    )

    equilibrium_points = np.array([point[1:] for point in equilibria])
    eigenvalue_text = []
    for _, x_value, y_value, z_value in equilibria:
        jacobian = np.array(
            [
                [-sigma, sigma, 0.0],
                [rho - z_value, -1.0, -x_value],
                [y_value, x_value, -beta],
            ]
        )
        eigenvalues = np.linalg.eigvals(jacobian)
        eigenvalue_text.append(
            ", ".join(
                f"{value.real:.3f}{value.imag:+.3f}i" for value in eigenvalues
            )
        )

    figure.add_trace(
        go.Scatter3d(
            x=equilibrium_points[:, 0],
            y=equilibrium_points[:, 1],
            z=equilibrium_points[:, 2],
            mode="markers+text",
            name="Equilibria",
            text=[point[0] for point in equilibria],
            textposition="top center",
            marker={"size": 5, "color": "#ffffff", "symbol": "diamond"},
            customdata=eigenvalue_text,
            hovertemplate=(
                "%{text}<br>x = %{x:.4f}<br>y = %{y:.4f}<br>z = %{z:.4f}"
                "<br>Jacobian eigenvalues: %{customdata}"
                "<extra>Equilibrium</extra>"
            ),
        )
    )

    figure.update_layout(
        template="plotly_dark",
        paper_bgcolor="#10151d",
        plot_bgcolor="#10151d",
        margin={"l": 0, "r": 0, "t": 24, "b": 0},
        height=720,
        legend={"x": 0.02, "y": 0.98},
        scene={
            "xaxis_title": "x",
            "yaxis_title": "y",
            "zaxis_title": "z",
            "aspectmode": "data",
            "camera": {"eye": {"x": 1.65, "y": 1.65, "z": 1.15}},
        },
    )

    output_path = Path(__file__).with_name("lorenz_attractor.html")
    plot_fragment = figure.to_html(
        full_html=False,
        include_plotlyjs=True,
        config={"responsive": True, "displaylogo": False},
    )

    template_path = Path(__file__).with_name("lorenz_template.html")
    page = template_path.read_text(encoding="utf-8")
    replacements = {
        "__PLOT_FRAGMENT__": plot_fragment,
        "__SIGMA__": f"{sigma:g}",
        "__RHO__": f"{rho:g}",
        "__BETA__": f"{beta:.4g}",
        "__INITIAL_STATE__": str(tuple(initial_state)),
        "__DT__": f"{dt:g}",
        "__FINAL_TIME__": f"{time[-1]:g}",
        "__CURVE_INTEGRAL__": curve_integral_latex,
    }
    for placeholder, value in replacements.items():
        page = page.replace(placeholder, value)

    output_path.write_text(page, encoding="utf-8")
    import webbrowser

    webbrowser.open(output_path.as_uri())
    print(f"Interactive graph saved to: {output_path}")


def main():
    dt = 0.005

    initial_state = (1.0, 1.0, 1.0)

    # Starting from around (1, 1, 1) gives the usual butterfly shape after a
    # short bit of settling in. The first few seconds change a lot depending on
    # the starting point, even if the long-term orbit ends up on the same attractor.
    trajectory = integrate_lorenz(initial_state=initial_state, dt=dt)

    # The two wings are built around unstable equilibrium points. The Jacobian
    # tells us what nearby points do locally, but the orbit keeps getting stretched
    # and folded, so that local linear picture does not stay accurate forever.
    # This is really a trajectory in phase space, not an actual butterfly moving
    # around in normal physical space.
    plot_attractor(trajectory, dt, initial_state)


if __name__ == "__main__":
    main()
