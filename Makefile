# Disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# Project name
EXECUTABLE := diablo

# Configuration settings
FC := h5pfc
FFLAGS := -fbacktrace -fbounds-check
LDFLAGS := -L$(HOME)/2decomp_fft/lib -l2decomp_fft -lfftw3
INCLUDE = -I$(HOME)/2decomp_fft/include
SRCDIR = src
OBJDIR = obj
MODDIR = mod
RM := rm -f

FC += $(FFLAGS)
FC += -J $(MODDIR)

# List of all source files
SRCS := $(wildcard src/*.f90)

# Create lists of the build artefacts in this project
OBJS := $(subst .f90,.o,$(SRCS))
OBJS := $(subst src,obj,$(OBJS))

# Declare all public targets
.PHONY: all clean
all: $(EXECUTABLE)

# Link the test executables
$(EXECUTABLE): $(OBJS) $(MODFILES)
# $(EXECUTABLE): diablo.o $(MODFILES)
	$(FC) $(LDFLAGS) -o $@ $^

# Create object files from Fortran source
$(OBJS): $(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(INCLUDE) -c -o $@ $<
# %.o: src/%.f90
# 	$(FC) $(INCLUDE) $(FFLAGS) $< -o $@

# Define all module interdependencies
obj/param.o: obj/grid.o
obj/diablo_io.o: obj/grid.o obj/param.o obj/fft.o obj/user_IC.o obj/periodic.o
obj/fft.o: obj/grid.o obj/param.o
obj/periodic.o: obj/param.o obj/fft.o obj/user_rhs.o
obj/diablo.o: obj/periodic.o obj/param.o obj/diablo_io.o obj/hdf5_mod.o
obj/hdf5_mod.o: obj/grid.o obj/param.o
obj/user_IC.o: obj/grid.o obj/param.o obj/periodic.o obj/fft.o
obj/user_rhs.o: obj/grid.o obj/param.o

clean:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod diablo