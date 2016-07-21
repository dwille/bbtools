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
LDINCS = -I $(CGNS_DIR)/include 
LDLIBS = -lm -L $(HDF5_DIR)/lib -L $(CGNS_DIR)/lib -lcgns -lhdf5

SRCC =	main.c 		\
	cgns_reader.c

EXTRA = Makefile	\
	main.h		\
	cgns_reader.h

# compile normally:
all: COPT += -O2
all: track-wave

# compile with debug flags
debug: COPT += -DDEBUG -g
debug: track-wave

# compile for batch job submission
batch: COPT += -O2 -DBATCH
batch: track-wave

# compile with profiler
profiler: COPT += -pg
profiler: track-wave

# compile with debug and batch job sumission
batchdebug: COPT += -DDEBUG -DBATCH -g
batchdebug: track-wave

OBJS = $(addprefix $(SRC_DIR)/, $(addsuffix .o, $(basename $(SRCC))))

$(OBJS):$(SRC_DIR)/%.o:$(SRC_DIR)/%.c
	$(GCC) $(COPT) -c $< $(LDINCS) -o $@

track-wave: $(OBJS)
	$(GCC) $(COPT) $+ $(LDLIBS) -o $(SIM_DIR)/track-wave

clean: 
	rm -f $(SRC_DIR)/*.o $(SIM_DIR)/track-wave
