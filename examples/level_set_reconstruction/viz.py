import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# User settings
# ============================================================

DATA_FILE = "build/level_set_viz/errors_master_Ellipsoid.txt"
OUTPUT_DIR = "build/level_set_viz/ellipsoid_plots"

# Radius used for convergence plots
CONVERGENCE_RADIUS = 2.5

# Expected convergence orders to show
REFERENCE_ORDERS = [1, 2]

# Save figures?
SAVE_FIGURES = True

# Show figures interactively?
SHOW_FIGURES = False

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ============================================================
# Load data
# ============================================================

df = pd.read_csv(DATA_FILE, sep=r"\s+")
# ============================================================
# Reconstruct weight column from sweep ordering
# ============================================================

weights = [
    "Wu4",
    "Wendland2",
    "Wendland4",
    "Wendland6",
    "Wu2",
]

runs_per_weight = 4 * 2 * 1 * 41   # Nx * methods * shapes * radii = 328

expected_rows = len(weights) * runs_per_weight

if len(df) != expected_rows:
    raise RuntimeError(
        f"Expected {expected_rows} rows, but found {len(df)}. "
        "Cannot safely reconstruct weight labels."
    )

df["weight"] = np.repeat(weights, runs_per_weight)

print(df[["nx", "method", "radius_cells", "weight"]].head())
print(df[["nx", "method", "radius_cells", "weight"]].tail())

print("Columns:")
print(df.columns.tolist())

print("\nNumber of cases:", len(df))

print("\nWeights:")
print(df["weight"].unique())

print("\nNx:")
print(np.sort(df["nx"].unique()))

print("\nMethods:")
print(df["method"].unique())

print("\nRadii:")
print(np.sort(df["radius_cells"].unique()))


# ============================================================
# Plot settings
# ============================================================

error_types = {
    "position_error": {
        "label": "Position Error",
        "orders": [2, 3],
    },
    "normal_error": {
        "label": "Normal Error",
        "orders": [1, 2],
    },
    "curvature_error": {
        "label": "Curvature Error",
        "orders": [0, 1],
    },
}

weights = [
    "Wu2",
    "Wu4",
    "Wendland2",
    "Wendland4",
    "Wendland6",
]

methods = ["Jibben", "LVIRA"]

nx_values = sorted(df["nx"].unique())

# ============================================================
# Visual styles
# ============================================================

nx_colors = {
    16:  "tab:blue",
    32:  "tab:orange",
    64:  "tab:green",
    128: "tab:red",
}

method_markers = {
    "Jibben": "o",
    "LVIRA":  "s",
}

method_linestyles = {
    "Jibben": "-",
    "LVIRA":  "--",
}
# ============================================================
# Plot 1:
# Error vs kernel radius
#
# Color  -> Nx
# Marker -> reconstruction method
# ============================================================

radius_output = os.path.join(OUTPUT_DIR, "radius_sweep")
os.makedirs(radius_output, exist_ok=True)

for weight in weights:

    weight_data = df[df["weight"] == weight]

    for error_column, error_info in error_types.items():

        error_label = error_info["label"]

        plt.figure(figsize=(8, 6))

        for nx in nx_values:

            for method in methods:

                subset = weight_data[
                    (weight_data["method"] == method)
                    & (weight_data["nx"] == nx)
                ].sort_values("radius_cells")

                if subset.empty:
                    continue

                plt.semilogy(
                    subset["radius_cells"],
                    subset[error_column],
                    color=nx_colors[nx],
                    marker=method_markers[method],
                    linestyle=method_linestyles[method],
                    linewidth=2.5,
                    markersize=5,
                    markevery=2,
                    label=rf"$N_x={nx}$, {method}",
                )

        plt.xlabel(r"Kernel Radius, $\delta/\Delta x$")
        plt.ylabel(error_label)

        plt.title(
            f"{error_label} vs Kernel Radius\n"
            f"{weight}"
        )

        plt.grid(True, which="both", alpha=0.3)
        plt.legend(
            fontsize=8,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
        )
        plt.tight_layout()

        if SAVE_FIGURES:
            filename = (
                f"{radius_output}/"
                f"{weight}_{error_column}_vs_radius.png"
            )

            plt.savefig(
                filename,
                dpi=300,
                bbox_inches="tight",
            )

        if SHOW_FIGURES:
            plt.show()

        plt.close()


# ============================================================
# Plot 2:
# Error convergence vs Nx
#
# One figure for each error type = 3 figures
#
# Color     -> weighting function
# Linestyle -> reconstruction method
#
# Data is taken at CONVERGENCE_RADIUS.
# ============================================================

convergence_output = os.path.join(OUTPUT_DIR, "convergence")
os.makedirs(convergence_output, exist_ok=True)

# ------------------------------------------------------------
# Visual styles
# ------------------------------------------------------------

weight_colors = {
    "Wu2":       "tab:blue",
    "Wu4":       "tab:orange",
    "Wendland2": "tab:green",
    "Wendland4": "tab:red",
    "Wendland6": "tab:purple",
}

method_linestyles = {
    "Jibben": "-",
    "LVIRA":  "--",
}

# ------------------------------------------------------------
# Select desired radius
# ------------------------------------------------------------

radius_data = df[
    np.isclose(
        df["radius_cells"],
        CONVERGENCE_RADIUS,
        rtol=0.0,
        atol=1.0e-10,
    )
]

if radius_data.empty:
    raise RuntimeError(
        f"No data found for radius = {CONVERGENCE_RADIUS}"
    )


# ============================================================
# Generate one convergence plot per error type
# ============================================================

for error_column, error_info in error_types.items():

    error_label = error_info["label"]
    reference_orders = error_info["orders"]

    plt.figure(figsize=(8, 6))

    plotted_data = []

    # --------------------------------------------------------
    # Numerical results
    # --------------------------------------------------------

    for weight in weights:

        for method in methods:

            subset = radius_data[
                (radius_data["weight"] == weight)
                & (radius_data["method"] == method)
            ].sort_values("nx")

            if subset.empty:
                continue

            x = subset["nx"].to_numpy()
            y = subset[error_column].to_numpy()

            plt.loglog(
                x,
                y,
                color=weight_colors[weight],
                linestyle=method_linestyles[method],
                marker="o",
                linewidth=2,
                markersize=6,
                label=f"{weight}, {method}",
            )

            plotted_data.append((x, y))

   # --------------------------------------------------------
    # Reference convergence curves
    #
    # First reference order  -> anchored to LVIRA
    # Second reference order -> anchored to Jibben
    # --------------------------------------------------------

    REFERENCE_WEIGHT = "Wu2"

    lvira_ref = radius_data[
        (radius_data["weight"] == REFERENCE_WEIGHT)
        & (radius_data["method"] == "LVIRA")
    ].sort_values("nx")

    jibben_ref = radius_data[
        (radius_data["weight"] == REFERENCE_WEIGHT)
        & (radius_data["method"] == "Jibben")
    ].sort_values("nx")


    # --------------------------------------------------------
    # LVIRA reference
    # --------------------------------------------------------

    lvira_order = reference_orders[0]

    if not lvira_ref.empty:

        x_reference = lvira_ref["nx"].to_numpy()
        y_data = lvira_ref[error_column].to_numpy()

        x0 = x_reference[0]
        y0 = y_data[0]

        y_reference = y0 * (x_reference / x0) ** (-lvira_order)

        plt.loglog(
            x_reference,
            y_reference,
            color="black",
            linestyle=":",
            linewidth=1.5,
            label=rf"$O(N_x^{{-{lvira_order}}})$",
        )


    # --------------------------------------------------------
    # Jibben reference
    # --------------------------------------------------------

    jibben_order = reference_orders[1]

    if not jibben_ref.empty:

        x_reference = jibben_ref["nx"].to_numpy()
        y_data = jibben_ref[error_column].to_numpy()

        x0 = x_reference[0]
        y0 = y_data[0]

        y_reference = y0 * (x_reference / x0) ** (-jibben_order)

        plt.loglog(
            x_reference,
            y_reference,
            color="black",
            linestyle="-.",
            linewidth=1.5,
            label=rf"$O(N_x^{{-{jibben_order}}})$",
        )

    # --------------------------------------------------------
    # Figure formatting
    # --------------------------------------------------------

    plt.xlabel(r"$N_x$")
    plt.ylabel(error_label)

    plt.title(
        f"{error_label} Convergence\n"
        rf"$\delta/\Delta x={CONVERGENCE_RADIUS}$"
    )

    plt.grid(True, which="both", alpha=0.3)

    # Legend outside plot on right
    plt.legend(
        fontsize=8,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        handlelength=4.0,
    )
    plt.tight_layout()

    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    if SAVE_FIGURES:

        filename = (
            f"{convergence_output}/"
            f"{error_column}_R{CONVERGENCE_RADIUS}.png"
        )

        plt.savefig(
            filename,
            dpi=300,
            bbox_inches="tight",
        )

    if SHOW_FIGURES:
        plt.show()

    plt.close()



# ============================================================
# Plot 3:
# Convergence vs Nx for different kernel radii
#
# One figure for each:
#     weight x method x error type
#
# 5 weights x 2 methods x 3 errors = 30 figures
#
# Each curve -> different kernel radius
# ============================================================

radius_convergence_output = os.path.join(
    OUTPUT_DIR, "radius_convergence"
)
os.makedirs(radius_convergence_output, exist_ok=True)


# ------------------------------------------------------------
# Radii to plot
#
# Choose whichever radii you want displayed.
# ------------------------------------------------------------

PLOT_RADII = np.arange(1.0, 5.01, 0.5)

# To instead plot every available radius, use:
# PLOT_RADII = np.sort(df["radius_cells"].unique())

# ------------------------------------------------------------
# Radius color mapping
#
# Blue -> small radius
# Red  -> large radius
# ------------------------------------------------------------

cmap = plt.get_cmap("coolwarm")

radius_min = np.min(PLOT_RADII)
radius_max = np.max(PLOT_RADII)

norm = plt.Normalize(
    vmin=radius_min,
    vmax=radius_max
)


# ------------------------------------------------------------
# Generate plots
# ------------------------------------------------------------

for weight in weights:

    for method in methods:

        case_data = df[
            (df["weight"] == weight)
            & (df["method"] == method)
        ]

        for error_column, error_info in error_types.items():

            error_label = error_info["label"]

            plt.figure(figsize=(8, 6))

            # ------------------------------------------------
            # One curve for each radius
            # ------------------------------------------------

            for radius in PLOT_RADII:

                subset = case_data[
                    np.isclose(
                        case_data["radius_cells"],
                        radius,
                        rtol=0.0,
                        atol=1.0e-10,
                    )
                ].sort_values("nx")

                if subset.empty:
                    continue

                x = subset["nx"].to_numpy()
                y = subset[error_column].to_numpy()

                # Radius determines color
                color = cmap(norm(radius))

                plt.loglog(
                    x,
                    y,
                    color=color,
                    marker="o",
                    linewidth=1.5,
                    markersize=5,
                    label=rf"$\delta/\Delta x={radius:.1f}$",
                )

            # ------------------------------------------------
            # Figure formatting
            # ------------------------------------------------

            plt.xlabel(r"$N_x$")
            plt.ylabel(error_label)

            plt.title(
                f"{error_label} Convergence\n"
                f"{weight}, {method}"
            )

            plt.grid(
                True,
                which="both",
                alpha=0.3,
            )

            # Legend outside plot
            sm = plt.cm.ScalarMappable(
                cmap=cmap,
                norm=norm
            )

            sm.set_array([])

            cbar = plt.colorbar(
                sm,
                ax=plt.gca(),
                pad=0.02
            )

            cbar.set_label(
                r"Kernel Radius, $\delta/\Delta x$"
            )
            plt.tight_layout()

            # ------------------------------------------------
            # Save
            # ------------------------------------------------

            if SAVE_FIGURES:

                filename = (
                    f"{radius_convergence_output}/"
                    f"{weight}_{method}_{error_column}.png"
                )

                plt.savefig(
                    filename,
                    dpi=300,
                    bbox_inches="tight",
                )

            if SHOW_FIGURES:
                plt.show()

            plt.close()

# ============================================================
# Plot 4:
# Optimal kernel radius vs Nx
#
# One figure for each error type = 3 figures
#
# x-axis    -> Nx (log scale)
# y-axis    -> radius producing minimum error
#
# Color     -> weighting function
# Linestyle -> reconstruction method
#
# Solid  -> Jibben
# Dashed -> LVIRA
# ============================================================

optimal_radius_output = os.path.join(
    OUTPUT_DIR, "optimal_radius"
)
os.makedirs(optimal_radius_output, exist_ok=True)


# ------------------------------------------------------------
# Visual styles
# ------------------------------------------------------------

weight_colors = {
    "Wu2":       "tab:blue",
    "Wu4":       "tab:orange",
    "Wendland2": "tab:green",
    "Wendland4": "tab:red",
    "Wendland6": "tab:purple",
}

method_linestyles = {
    "Jibben": "-",
    "LVIRA": "--",
}


# ============================================================
# Generate one optimal-radius plot per error type
# ============================================================

for error_column, error_info in error_types.items():

    error_label = error_info["label"]

    plt.figure(figsize=(8, 6))

    # --------------------------------------------------------
    # Loop over weighting functions and methods
    # --------------------------------------------------------

    for weight in weights:

        for method in methods:

            case_data = df[
                (df["weight"] == weight)
                & (df["method"] == method)
            ]

            nx_values = []
            optimal_radii = []

            # ------------------------------------------------
            # Find optimal radius independently at each Nx
            # ------------------------------------------------

            for nx in sorted(case_data["nx"].unique()):

                subset = case_data[
                    case_data["nx"] == nx
                ]

                if subset.empty:
                    continue

                # Index of row having minimum error
                min_index = subset[error_column].idxmin()

                # Corresponding optimal radius
                optimal_radius = subset.loc[
                    min_index, "radius_cells"
                ]

                nx_values.append(nx)
                optimal_radii.append(optimal_radius)

            # ------------------------------------------------
            # Plot optimal radius vs Nx
            # ------------------------------------------------

            plt.semilogx(
                nx_values,
                optimal_radii,
                color=weight_colors[weight],
                linestyle=method_linestyles[method],
                marker="o",
                linewidth=2,
                markersize=6,
                label=f"{weight}, {method}",
            )


    # --------------------------------------------------------
    # Figure formatting
    # --------------------------------------------------------

    plt.xlabel(r"$N_x$")
    plt.ylabel(
        rf"Optimal Kernel Radius, $\delta/\Delta x$"
    )

    plt.title(
        f"Optimal Kernel Radius Based on {error_label}"
    )

    plt.grid(
        True,
        which="both",
        alpha=0.3,
    )

    # Since Nx values are specifically 16, 32, 64, 128
    # plt.xticks(
    #     [16, 32, 64, 128],
    #     ["16", "32", "64", "128"],
    # )

    # The sweep covers radius = 1.0 -> 5.0
    plt.ylim(0.9, 5.1)

    # --------------------------------------------------------
    # Legend
    # --------------------------------------------------------

    plt.legend(
        fontsize=8,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        handlelength=4.0,
    )

    plt.tight_layout()


    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    if SAVE_FIGURES:

        filename = (
            f"{optimal_radius_output}/"
            f"{error_column}_optimal_radius.png"
        )

        plt.savefig(
            filename,
            dpi=300,
            bbox_inches="tight",
        )

    if SHOW_FIGURES:
        plt.show()

    plt.close()


print("Plot 4 complete.")


print("\nPlotting complete.")
print(f"Figures saved to: {OUTPUT_DIR}")