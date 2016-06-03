## DEPENDENCIES ##
HDF5_DIR = /home-1/dwillen3@jhu.edu/.deps/hdf5
CGNS_DIR = /home-1/dwillen3@jhu.edu/.deps/cgns

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
all: partVel

# compile with debug flags
debug: COPT += -DDEBUG -g
debug: partVel

# compile for batch job submission
batch: COPT += -O2 -DBATCH
batch: partVel

# compile with profiler
profiler: COPT += -pg
profiler: partVel

OBJS = $(addprefix $(SRC_DIR)/, $(addsuffix .o, $(basename $(SRCC))))

$(OBJS):$(SRC_DIR)/%.o:$(SRC_DIR)/%.c
	$(GCC) $(COPT) -c $< $(LDINCS) -o $@

partVel: $(OBJS)
	$(GCC) $(COPT) $+ $(LDLIBS) -o $(SIM_DIR)/partVel

clean: 
	rm -f $(SRC_DIR)/*.o $(SIM_DIR)/partVel
