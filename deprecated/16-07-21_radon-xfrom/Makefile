## DEPENDENCIES ##
HDF5_DIR = /opt/hdf5
CGNS_DIR = /opt/cgns

## COMPILERS
GCC = gcc

SIM_DIR = sim
SRC_DIR = src

COPT = -std=c99 -pedantic -Wall -Wextra
LDINCS = -I $(CGNS_DIR)/include
LDLIBS = -lm -L $(HDF5_DIR)/lib -L $(CGNS_DIR)/lib -lcgns -lhdf5

SRCC =	main.c 		\
	radon.c		\
	cgns_reader.c

EXTRA = Makefile	\
	main.h		\
	radon.h		\
	cgns_reader.h

# compile normally:
all: COPT += -O2
all: radon-xform

# compile with debug flags
debug: COPT += -DDEBUG -g
debug: radon-xform

# compile with profiler
profiler: COPT += -pg
profiler: radon-xform

OBJS = $(addprefix $(SRC_DIR)/, $(addsuffix .o, $(basename $(SRCC))))

$(OBJS):$(SRC_DIR)/%.o:$(SRC_DIR)/%.c
	$(GCC) $(COPT) -c $< $(LDINCS) -o $@

radon-xform: $(OBJS)
	$(GCC) $(COPT) $+ $(LDLIBS) -o $(SIM_DIR)/radon-xform

clean: 
	rm -f $(SRC_DIR)/*.o $(SIM_DIR)/radon-xform
