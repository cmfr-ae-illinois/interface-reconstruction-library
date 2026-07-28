from pathlib import Path

# params to change -------------------------------------------
# ------------------------------------------------------------
# case_name = "deformation3d"
case_name = "rotation3d"
reconstruction_names = [
    "lvira",
    "vf",
    "vf2",
    "cf",
    "plicnet",
]
max_levels = [
    3,
    4,
    5,
    6,
    7,
]
cfl = 0.25
output_dir = Path(
    "/home/parinht2/Repositories/"
    "interface-reconstruction-library/"
    "examples/amrex_advector/input_files"
)
checkpoint_path = "."
interface_output_path = "."
# ------------------------------------------------------------
# ------------------------------------------------------------

max_step = 100000
advection_name = "remap"
is_periodic = [1, 1, 1]
amr_v = 1 # verbosity
n_cell = [2, 2, 2] # base resolution
ref_ratio_value = 2 # Refinement ratio between every level
blocking_factor_x = 2
blocking_factor_y = 2
blocking_factor_z = 2
max_grid_size = 8
regrid_int = 1
do_reflux = 1
# plot_file = "plt"
plot_int = -1
chk_int = -1

# CASE-SPECIFIC SETTINGS
def get_case_settings(case):
    if case == "deformation3d":
        stop_time = 3.0
        prob_lo = [0.0, 0.0, 0.0,]
        prob_hi = [1.0, 1.0, 1.0,]
    elif case == "rotation3d":
        stop_time = 1.0
        prob_lo = [-0.5, -0.5, -0.5,]
        prob_hi = [0.5, 0.5, 0.5,]
    else:
        raise ValueError(
            f"Unknown case_name: '{case}'\n"
            "Supported cases are:\n"
            "  deformation3d\n"
            "  rotation3d"
        )
    return stop_time, prob_lo, prob_hi



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
    return str(value)


# GENERATE INPUT FILE
def generate_input_file(
    reconstruction_name,
    max_level,
    stop_time,
    prob_lo,
    prob_hi,
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
amr.plot_int   = {plot_int}    # number of timesteps between plot files
                        # if negative then no plot files will be written
amr.interface_output_path = {interface_output_path}

# *****************************************************************
# Checkpoint name and frequency
# *****************************************************************
amr.chk_file = {chk_filename}      # root name of checkpoint file
amr.chk_int  = {chk_int}       # number of timesteps between checkpoint files
                        # if negative then no checkpoint files will be written
amr.checkpoint_path = {checkpoint_path}
"""

    # Write file
    input_path.write_text(content)
    print(f"Generated: {input_path}")
    print(f"  case             = {case_name}")
    print(f"  reconstruction   = {reconstruction_name}")
    print(f"  CFL              = {cfl}")
    print(f"  max_level        = {max_level}")
    print(f"  ref_ratio        = {format_list(ref_ratios)}")
    print(f"  stop_time        = {stop_time}")
    print(f"  prob_lo          = {format_list(prob_lo)}")
    print(f"  prob_hi          = {format_list(prob_hi)}")
    print(f"  checkpoint root  = {chk_filename}")
    print(f"  checkpoint path  = {checkpoint_path}")
    print(f"  interface path   = {interface_output_path}")
    print()

def main():
    stop_time, prob_lo, prob_hi = get_case_settings(case_name)
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
            )
            num_files += 1
# RUN
if __name__ == "__main__":
    main()