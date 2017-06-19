## DEPENDENCIES ##
ifeq ($(USER),dwille)
  $(info User is $(USER))
  #HDF5_DIR = /opt/hdf5
  #CGNS_DIR = /opt/cgns
  HDF5_DIR = /home/dwille/.deps/phdf5/1.8.18
  CGNS_DIR = /home/dwille/.deps/pcgns/3.3.0
else ifeq ($(USER),dwillen3@jhu.edu)
  $(info USER is dwillen3@jhu.edu)
  HDF5_DIR = /home-1/dwillen3@jhu.edu/.deps/phdf5/1.8.18
  CGNS_DIR = /home-1/dwillen3@jhu.edu/.deps/pcgns/3.3.0
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
all: f-rec-1D-part

# compile with debug flags
debug: COPT += -DDEBUG -g
debug: f-rec-1D-part

# compile with profiler
profiler: COPT += -pg
profiler: f-rec-1D-part

OBJS = $(addprefix $(SRC_DIR)/, $(addsuffix .o, $(basename $(SRCC))))

$(OBJS):$(SRC_DIR)/%.o:$(SRC_DIR)/%.c
	$(GCC) $(COPT) -c $< $(LDINCS) -o $@

f-rec-1D-part: $(OBJS)
	$(GCC) $(COPT) $+ $(LDLIBS) -o $(SIM_DIR)/f-rec-1D-part

clean: 
	rm -f $(SRC_DIR)/*.o $(SIM_DIR)/f-rec-1D-part
