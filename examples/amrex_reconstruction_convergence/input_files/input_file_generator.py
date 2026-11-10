import argparse
from pathlib import Path

shape = "orthocircle"
reconstruction_methods = [
    "lvira","vf",
]

binary_file = Path(
    "/home/parinht2/Repositories/interface-reconstruction-library/build/"
    "reconstruction_binary_files/orthocircle_no_s2/"
    "orthocircle_moments_Nx1024.bin"
)
nx_fine = 1024
factors = [4, 8, 16, 32, 64]
max_grid_size = 32

volume_moment_order = 2
surface_moment_order = 1
do_volume = 1
do_surface = 1
do_curvature = 1

output_interface = 1
output_csv = 0
output_cell_errors = 0

output_root_dir = Path(
    "/home/parinht2/Repositories/interface-reconstruction-library/build/"
    "reconstruction_convergence"
)

input_filename = "inputs_{shape}_{method}"
csv_filename = "{shape}_{method}_metrics.csv"
cell_error_output_name = "{shape}_{method}_cell_errors"
# -----------------------------------------------------------------------------


def format_list(values):
    """Return a sequence in AMReX's space-separated input format."""
    return " ".join(str(value) for value in values)


def expand_setting(value, method):
    """Expand {shape} and {method} placeholders in a string or Path setting."""
    return str(value).format(shape=shape, method=method)


def validate_filename_part(value, setting_name):
    """Reject values that could unexpectedly create nested output paths."""
    if not value or Path(value).name != value or value in {".", ".."}:
        raise ValueError(
            f"{setting_name} must be a non-empty filename component; got {value!r}"
        )


def generate_input_file(method, output_dir):
    """Generate one input file for a reconstruction method."""
    validate_filename_part(shape, "shape")
    validate_filename_part(method, "reconstruction method")

    current_input_filename = expand_setting(input_filename, method)
    current_csv_filename = expand_setting(csv_filename, method)
    cell_error_name = expand_setting(cell_error_output_name, method)
    validate_filename_part(current_input_filename, "input filename")
    validate_filename_part(current_csv_filename, "CSV filename")
    validate_filename_part(cell_error_name, "cell-error output name")

    input_path = output_dir / current_input_filename
    csv_file = output_dir / current_csv_filename
    content = f"""# implicit surface
shape = {shape}

# reconstruction method
reconstruction_method = {method}

# binary file containing the true moments at the finest resolution
binary_file = \"{binary_file}\"

# base grid size
nx_fine = {nx_fine}
# refinement factors
factors = {format_list(factors)}

max_grid_size = {max_grid_size}

# moment orders for reconstruction
volume_moment_order = {volume_moment_order}
surface_moment_order = {surface_moment_order}
do_volume = {do_volume}
do_surface = {do_surface}

# mean curvature error at the cell center projected onto the exact surface
do_curvature = {do_curvature}

# to output reconstructed interface
output_interface = {output_interface}
interface_output_dir = \"{output_dir}\"

# to output csv with enabled volumetric, surface, and curvature error norms
output_csv = {output_csv}
csv_file = \"{csv_file}\"

# to output cell errors for volumetric and surface moments
output_cell_errors = {output_cell_errors}
cell_error_output_dir = \"{output_dir}\"
cell_error_output_name = {cell_error_name}
"""

    input_path.write_text(content, encoding="utf-8")
    print(f"Generated: {input_path}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate one reconstruction-convergence input per method."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=output_root_dir,
        help=(
            "root directory for generated inputs and run outputs "
            f"(default: {output_root_dir})"
        ),
    )
    parser.add_argument(
        "--reconstruction-methods",
        nargs="+",
        metavar="METHOD",
        default=reconstruction_methods,
        help="methods to generate (default: the reconstruction_methods list)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir = args.output_dir.expanduser().resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for method in args.reconstruction_methods:
        generate_input_file(method, args.output_dir)


if __name__ == "__main__":
    main()
