## DEPENDENCIES ##
ifeq ($(USER),dwille)
  $(info User is $(USER))
  HDF5_DIR = /opt/hdf5
  CGNS_DIR = /opt/cgns
else ifeq ($(USER),dwillen3@jhu.edu)
  $(info USER is dwillen3@jhu.edu)
  HDF5_DIR = /home-1/dwillen3@jhu.edu/.deps/hdf5
  CGNS_DIR = /home-1/dwillen3@jhu.edu/.deps/cgns
  FFTW_DIR = /home-1/dwillen3@jhu.edu/.deps/fftw-3.3.4
endif	

## COMPILERS
GCC = gcc

SIM_DIR = sim
SRC_DIR = src

COPT = -std=c99 -pedantic -Wall -Wextra

ifeq ($(USER),dwille)
  LDINCS = -I $(CGNS_DIR)/include 
  LDLIBS = -lfftw3 -lm -L $(HDF5_DIR)/lib -L $(CGNS_DIR)/lib -lcgns -lhdf5
else ifeq ($(USER),dwillen3@jhu.edu)
  LDINCS = -I $(FFTW_DIR)/include -I $(CGNS_DIR)/include 
  LDLIBS = -L $(FFTW_DIR)/lib -lfftw3 -lm -L $(HDF5_DIR)/lib -L $(CGNS_DIR)/lib -lcgns -lhdf5
endif	

SRCC =	main.c 		\
	cgns_reader.c

EXTRA = Makefile	\
	main.h		\
	cgns_reader.h

# compile normally:
all: COPT += -O2
all: f-rec-1D-flow

# compile with debug flags
debug: COPT += -DDEBUG -g
debug: f-rec-1D-flow

# compile with profiler
profiler: COPT += -pg
profiler: f-rec-1D-flow

OBJS = $(addprefix $(SRC_DIR)/, $(addsuffix .o, $(basename $(SRCC))))

$(OBJS):$(SRC_DIR)/%.o:$(SRC_DIR)/%.c
	$(GCC) $(COPT) -c $< $(LDINCS) -o $@

f-rec-1D-flow: $(OBJS)
	$(GCC) $(COPT) $+ $(LDLIBS) -o $(SIM_DIR)/f-rec-1D-flow #-debug/f-rec-1D-flow

clean: 
	rm -f $(SRC_DIR)/*.o $(SIM_DIR)/f-rec-1D-flow #-debug/f-rec-flow-1D
