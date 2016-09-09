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

ifeq ($(USER),dwille)
  LDINCS = -I $(CGNS_DIR)/include 
  LDLIBS = -L $(HDF5_DIR)/lib -L $(CGNS_DIR)/lib -lcgns -lhdf5
else ifeq ($(USER),dwillen3@jhu.edu)
  LDINCS = -I $(CGNS_DIR)/include 
  LDLIBS = -L $(HDF5_DIR)/lib -L $(CGNS_DIR)/lib -lcgns -lhdf5
endif	

SRCC =	main.c 		\
	cgns_reader.c

EXTRA = Makefile	\
	main.h		\
	cgns_reader.h

# compile normally:
all: COPT += -O2
all: part-hist

# compile with debug flags
debug: COPT += -DDEBUG -g
debug: part-hist

# compile with profiler
profiler: COPT += -pg
profiler: part-hist

OBJS = $(addprefix $(SRC_DIR)/, $(addsuffix .o, $(basename $(SRCC))))

$(OBJS):$(SRC_DIR)/%.o:$(SRC_DIR)/%.c
	$(GCC) $(COPT) -c $< $(LDINCS) -o $@

part-hist: $(OBJS)
	$(GCC) $(COPT) $+ $(LDLIBS) -o $(SIM_DIR)/part-hist

clean: 
	rm -f $(SRC_DIR)/*.o $(SIM_DIR)/part-hist
