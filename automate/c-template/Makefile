# Dependencies
HDF5_DIR = /home-1/dwillen3@jhu.edu/.deps/hdf5
CGNS_DIR = /home-1/dwillen3@jhu.edu/.deps/cgns
GCC = gcc

SRC_DIR = src
SIM_DIR = .

# make -g -O2 to optimize
COPT = -std=c99 -pedantic -Wall -Wextra
LDINCS = -I $(CGNS_DIR)/include
LDLIBS = -lm -L $(HDF5_DIR)/lib -L $(CGNS_DIR)/lib -lcgns -lhdf5

SRCC =	template.c 	\
	reader.c

EXTRA = Makefile	\
	template.h	\
	reader.h	\

# compile normally:
all: COPT += -O2
all: template

# compile with debug flags
debug: COPT += -DDEBUG -g
debug: template

OBJS = $(addprefix $(SRC_DIR)/, $(addsuffix .o, $(basename $(SRCC))))

$(OBJS):$(SRC_DIR)/%.o:$(SRC_DIR)/%.c
	$(GCC) $(COPT) -c $< $(LDINCS) -o $@

template: $(OBJS)
	$(GCC) $(COPT) $+ $(LDLIBS) -o $(SIM_DIR)/template

clean:
	rm $(SRC_DIR)/*.o $(SIM_DIR)/template
