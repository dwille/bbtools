## DEPENDENCIES ##
HDF5_DIR = /opt/hdf5
CGNS_DIR = /opt/cgns

## COMPILERS
GCC = gcc

SIM_DIR = sim
SRC_DIR = src

COPT = -std=c99 -pedantic -Wall -Wextra
LDINCS = -I $(CGNS_DIR)/include 
LDLIBS = -lfftw3 -lm -L $(HDF5_DIR)/lib -L $(CGNS_DIR)/lib -lcgns -lhdf5

SRCC =	main.c 		\
	cgns_reader.c

EXTRA = Makefile	\
	main.h		\
	cgns_reader.h

# compile normally:
all: COPT += -O2
all: f-rec-flow-1D

# compile with debug flags
debug: COPT += -DDEBUG -g
debug: f-rec-flow-1D

# compile for batch job submission
batch: COPT += -O2 -DBATCH
batch: f-rec-flow-1D

# compile with profiler
profiler: COPT += -pg
profiler: f-rec-flow-1D

# compile with debug and batch job sumission
batchdebug: COPT += -DDEBUG -DBATCH -g
batchdebug: f-rec-flow-1D

OBJS = $(addprefix $(SRC_DIR)/, $(addsuffix .o, $(basename $(SRCC))))

$(OBJS):$(SRC_DIR)/%.o:$(SRC_DIR)/%.c
	$(GCC) $(COPT) -c $< $(LDINCS) -o $@

f-rec-flow-1D: $(OBJS)
	$(GCC) $(COPT) $+ $(LDLIBS) -o $(SIM_DIR)/f-rec-flow-1D #-debug/f-rec-flow-1D

clean: 
	rm -f $(SRC_DIR)/*.o $(SIM_DIR)/f-rec-flow-1D #-debug/f-rec-flow-1D
