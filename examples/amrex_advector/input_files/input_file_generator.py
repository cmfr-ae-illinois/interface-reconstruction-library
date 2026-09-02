from pathlib import Path

# ------------------------------------------------------------
# params to change 
# ------------------------------------------------------------
case_name = "deformation3d"
# case_name = "rotation3d"
reconstruction_names = [
    "pu",
]
max_levels = [
    3,
]
cfl = 1.00
output_dir = Path(
    "/home/parinht2/Repositories/"
    "interface-reconstruction-library/"
    "examples/amrex_advector/input_files"
)
checkpoint_path = "/home/parinht2/Repositories/interface-reconstruction-library/build/temporary"
interface_output_path = "/home/parinht2/Repositories/interface-reconstruction-library/build/temporary"
plot_times = [0.0, 0.25, 0.5, 0.75, 1.0]
chk_times = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
# ------------------------------------------------------------
# ------------------------------------------------------------

max_step = 10000
advection_name = "remap"
is_periodic = [1, 1, 1]
amr_v = 1
n_cell = [8, 8, 8]
ref_ratio_value = 2
blocking_factor_x = 8
blocking_factor_y = 8
blocking_factor_z = 8
max_grid_size = 32
regrid_int = 1
do_reflux = 1
plot_int = -1
chk_int = -1

# CASE-SPECIFIC SETTINGS
def get_case_settings(case):
    if case == "deformation3d":
        stop_time = 3.0
        prob_lo = [0.0, 0.0, 0.0,]
        prob_hi = [1.0, 1.0, 1.0,]
        velocity_field = 1
    elif case == "rotation3d":
        stop_time = 1.0
        prob_lo = [-0.5, -0.5, -0.5,]
        prob_hi = [0.5, 0.5, 0.5,]
        velocity_field = 1
    else:
        raise ValueError(
            f"Unknown case_name: '{case}'\n"
            "Supported cases are:\n"
            "  deformation3d\n"
            "  rotation3d"
        )
    return stop_time, prob_lo, prob_hi, velocity_field



# FORMATTING FUNCTIONS
def format_list(values):
    """
    Python list to a space-separated AMReX input string.
    """
    return " ".join(str(value) for value in values)


def format_cfl_for_filename(value):
    """
    Format CFL for use in filenames.
    """
    return f"{value:.2f}"


def format_optional_input(name, values):
    if values is None:
        return f"# {name} = 0.0 1.0"
    return f"{name} = {format_list(values)}"


# GENERATE INPUT FILE
def generate_input_file(
    reconstruction_name,
    max_level,
    stop_time,
    prob_lo,
    prob_hi,
    velocity_field,
):
    ref_ratios = [ref_ratio_value] * max_level
    cfl_string = format_cfl_for_filename(cfl)
    input_case_identifier = (
            f"{case_name}_"
            f"{reconstruction_name}_"
            f"cfl{cfl_string}_"
            f"{max_level}"
        )
    case_identifier = (
        f"{case_name}_"
        f"{reconstruction_name}_"
        f"cfl{cfl_string}_"
        f"{max_level}_"
    )
    # Input filename
    input_filename = f"input_{input_case_identifier}"
    # Checkpoint root name written inside the input file
    chk_filename = f"chk_{case_identifier}"
    # plot file name
    plot_filename = f"plt_{case_identifier}"
    # interface PVD name follows plot/checkpoint naming without the trailing factor-style underscore
    interface_pvd_file = f"interface_{input_case_identifier}.pvd"
    # Full path
    input_path = output_dir / input_filename

    # Build input file contents
    content = f"""# *****************************************************************
# Run until nsteps == max_step or time == stop_time,
# whichever comes first
# *****************************************************************
max_step  = {max_step}
stop_time = {stop_time}

# *****************************************************************
# Which case to run?
# *****************************************************************
case.name = {case_name}

# *****************************************************************
# Which interface reconstruction method to use?
# *****************************************************************
reconstruction.name = {reconstruction_name}

# *****************************************************************
# Which interface advection method to use?
# *****************************************************************
advection.name = {advection_name}
adv.velocity_field = {velocity_field}  # 0 exact case velocity, 1 interpolated face velocity

# *****************************************************************
# Are we restarting from an existing checkpoint file?
# *****************************************************************
#amr.restart  = chk00060 # restart from this checkpoint file

# *****************************************************************
# Problem size and geometry
# *****************************************************************
geometry.prob_lo     = {format_list(prob_lo)}
geometry.prob_hi     = {format_list(prob_hi)}

# 1 = periodic
# 0 = walls (homogeneous Neumann, dphi/dn=0)
geometry.is_periodic = {format_list(is_periodic)}

# *****************************************************************
# VERBOSITY
# *****************************************************************
amr.v              = {amr_v}

# *****************************************************************
# Resolution and refinement
# *****************************************************************
amr.n_cell          = {format_list(n_cell)}
amr.max_level       = {max_level}       # maximum level number allowed --
                              # number of levels = max_level + 1

amr.ref_ratio       = {format_list(ref_ratios)} # refinement ratio between levels

# *****************************************************************
# Control of grid creation
# *****************************************************************
# Blocking factor for grid creation in each dimension --
# this ensures that every grid is coarsenable by a factor of 8 --
# this is mostly relevant for multigrid performance
amr.blocking_factor_x = {blocking_factor_x}
amr.blocking_factor_y = {blocking_factor_y}
amr.blocking_factor_z = {blocking_factor_z}

amr.max_grid_size   = {max_grid_size}

amr.regrid_int      = {regrid_int}       # how often to regrid

# *****************************************************************
# Time step control
# *****************************************************************
adv.cfl            = {cfl}     # CFL constraint for explicit advection

# *****************************************************************
# Should we reflux at coarse-fine boundaries?
# *****************************************************************
adv.do_reflux = {do_reflux}

# *****************************************************************
# Plotfile name and frequency
# *****************************************************************
amr.plot_file  = {plot_filename}    # root name of plot file
# amr.plot_dir = .
amr.plot_int   = {plot_int}    # number of timesteps between plot files
                        # if negative then no plot files will be written
{format_optional_input("amr.plot_times", plot_times)}
                        # fractions of stop_time for plot files
                        # ignored when amr.plot_int > 0
# amr.interface_output_path = {interface_output_path}
amr.interface_pvd_file = {interface_pvd_file}

# *****************************************************************
# Checkpoint name and frequency
# *****************************************************************
amr.chk_file = {chk_filename}      # root name of checkpoint file
# amr.chk_dir = .
amr.chk_int  = {chk_int}       # number of timesteps between checkpoint files
                        # if negative then no checkpoint files will be written
{format_optional_input("amr.chk_times", chk_times)}
                        # fractions of stop_time for checkpoint files
                        # ignored when amr.chk_int > 0
# amr.checkpoint_path = {checkpoint_path}
"""

    # Write file
    input_path.write_text(content)
    print(f"Generated: {input_path}")
    # print(f"  case             = {case_name}")
    # print(f"  reconstruction   = {reconstruction_name}")
    # print(f"  CFL              = {cfl}")
    # print(f"  max_level        = {max_level}")
    # print(f"  ref_ratio        = {format_list(ref_ratios)}")
    # print(f"  stop_time        = {stop_time}")
    # print(f"  prob_lo          = {format_list(prob_lo)}")
    # print(f"  prob_hi          = {format_list(prob_hi)}")
    # print(f"  velocity_field   = {velocity_field}")
    # print(f"  plot_times       = {format_list(plot_times) if plot_times else None}")
    # print(f"  chk_times        = {format_list(chk_times) if chk_times else None}")
    # print(f"  checkpoint root  = {chk_filename}")
    # print(f"  checkpoint path  = {checkpoint_path}")
    # print(f"  interface path   = {interface_output_path}")
    # print(f"  interface pvd    = {interface_pvd_file}")
    print()

def main():
    stop_time, prob_lo, prob_hi, velocity_field = get_case_settings(case_name)
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )
    num_files = 0
    for reconstruction_name in reconstruction_names:
        for max_level in max_levels:
            generate_input_file(
                reconstruction_name=reconstruction_name,
                max_level=max_level,
                stop_time=stop_time,
                prob_lo=prob_lo,
                prob_hi=prob_hi,
                velocity_field=velocity_field,
            )
            num_files += 1
# RUN
if __name__ == "__main__":
    main()
