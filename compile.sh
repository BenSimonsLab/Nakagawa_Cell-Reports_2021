# The simulation can be compiled using various preprocessor flags:
#   - "labelx": clonally labels all X cells after equilibration of the system
#   - "labely": clonally labels all Y cells after equilibration of the system
#   - "labelz": clonally labels all Z cells after equilibration of the system
#   - "labelz2": clonally labels a fraction of Z cells (specified in the source code)
#                after equilibration of the system
#   - "labelcomp": labels all cells with a compartment-specific label
#   - "single": outputs all time-dependent joint clone size distributions
#               (can be combined with other flags)

# compilation using gfortran using fastest optimisation:
gfortran -cpp -Dlabelx -Dsingle germline.f90 -ffree-line-length-none -fbackslash -Ofast -o germline
