# Dependencies
CGNS_DIR = /opt/cgns
GCC = gcc

SRC_DIR = src
SIM_DIR = .

# make -g -O2 to optimize
COPT = -std=c99 -pedantic -Wall -Wextra
LDINCS = -I $(CGNS_DIR)/include
LDLIBS = -lm -L $(CGNS_DIR)/lib -lcgns

SRCC =	multipart.c 	\
	reader.c

EXTRA = Makefile	\
	multipart.h	\
	reader.h	\

# compile normally:
all: COPT += -O2
all: multipart

# compile with debug flags
debug: COPT += -DDEBUG -g
debug: multipart

OBJS = $(addprefix $(SRC_DIR)/, $(addsuffix .o, $(basename $(SRCC))))

$(OBJS):$(SRC_DIR)/%.o:$(SRC_DIR)/%.c
	$(GCC) $(COPT) -c $< $(LDINCS) -o $@

multipart: $(OBJS)
	$(GCC) $(COPT) $+ $(LDLIBS) -o $(SIM_DIR)/multipart

clean:
	rm $(SRC_DIR)/*.o $(SIM_DIR)/multipart
