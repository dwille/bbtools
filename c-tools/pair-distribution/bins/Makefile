## DEPENDENCIES ##
ifeq ($(USER),dwille)
  $(info User is $(USER))
  HDF5_DIR = /opt/hdf5
  CGNS_DIR = /opt/cgns
else ifeq ($(USER),dwillen3@jhu.edu)
  $(info USER is dwillen3@jhu.edu)
  HDF5_DIR = /home-1/dwillen3@jhu.edu/.deps/hdf5
  CGNS_DIR = /home-1/dwillen3@jhu.edu/.deps/cgns
endif	

## COMPILERS
GCC = gcc

SIM_DIR = sim
SRC_DIR = src

COPT = -std=c99 -pedantic -Wall -Wextra
LDINCS = -I $(CGNS_DIR)/include
LDLIBS = -lm -L $(HDF5_DIR)/lib -L $(CGNS_DIR)/lib -lcgns -lhdf5

SRCC =	main.c 		\
	reconstruct.c	\
	cgns_reader.c

EXTRA = Makefile	\
	main.h		\
	reconstruct.h	\
	cgns_reader.h

# compile normally:
all: COPT += -O2
all: part-pair-bin

# compile with debug flags
debug: COPT += -DDEBUG -g
debug: part-pair-bin

# compile with profiler
profiler: COPT += -O2 -pg
profiler: part-pair-bin

# compile for batch job submission
batch: COPT += -O2 -DBATCH
batch: part-pair-bin

# compile with debug and batch job sumission
batchdebug: COPT += -DDEBUG -DBATCH -g
batchdebug: part-pair-bin

OBJS = $(addprefix $(SRC_DIR)/, $(addsuffix .o, $(basename $(SRCC))))

$(OBJS):$(SRC_DIR)/%.o:$(SRC_DIR)/%.c
	$(GCC) $(COPT) -c $< $(LDINCS) -o $@

part-pair-bin: $(OBJS)
	$(GCC) $(COPT) $+ $(LDLIBS) -o $(SIM_DIR)/part-pair-bin

clean: 
	rm -f $(SRC_DIR)/*.o $(SIM_DIR)/part-pair-bin
