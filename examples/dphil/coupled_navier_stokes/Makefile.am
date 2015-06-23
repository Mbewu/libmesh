example_name  = picard_navier_stokes_read_mesh_1
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/dphil/picard_navier_stokes_read_mesh_1
data          = lid_driven_cavity_boundary_function.h  \
								moving_pipe_boundary_function.h \
								inflow_boundary_function.h \
								exodusII_io_helper_extended.h \
								exodusII_io_helper_extended.C \
								exodusII_io_extended.h \
								exodusII_io_extended.C \
								ball_boundary_flagging.h \
								xyz_position_flagging.h \
								z_position_flagging.h \
								ns_assembler_3d.h \
								ns_assembler_3d.C \
								picard.h \
								picard.C \
								optimised_stabilised_assembler_3d.h \
								optimised_stabilised_assembler_3d.C \
								navier_stokes_assembler.h \
								navier_stokes_assembler.C \
								coupled_assembler.h \
								augment_sparsity_on_interface.h \
								augment_sparsity_on_interface.C \
								surface_boundary.h \
								airway.h \
								airway.C \
								particle.h \
								particle.C \
								hofmann_particle_deposition.h \
								hofmann_particle_deposition.C \
								coupled_navier_stokes.h \
								coupled_navier_stokes.C \
								coupled_navier_stokes_1d_functions.C \
								coupled_navier_stokes_3d_functions.C \
								run.sh
sources       = $(data) run.sh

CLEANFILES = out_*.e

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
