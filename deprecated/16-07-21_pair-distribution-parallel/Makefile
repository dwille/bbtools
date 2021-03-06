## DEPENDENCIES ##
HDF5_DIR = /opt/hdf5
CGNS_DIR = /opt/cgns
CUSP_DIR = /opt/cusplibrary
CUDA_DIR = /opt/cuda

## COMPILERS
GCC = gcc
NVCC = nvcc

SIM_DIR = sim
SRC_DIR = src

COPT = -std=c99 -pedantic -Wall -Wextra
LDINCS = -I $(CGNS_DIR)/include -I $(CUSP_DIR)
LDLIBS = -lm -L $(HDF5_DIR)/lib -L $(CGNS_DIR)/lib -lcgns -lhdf5

CUDAOPT = -arch=sm_30 -Xcompiler -m64

CUDALIBS = -L $(CUDA_DIR)/lib64 -lcudart

SRCC =	main.c 		\
	reader.c

SRCCUDA = my_cuda.cu		\
	  cuda_kernel.cu

EXTRA = Makefile		\
	main.h		\
	cuda_kernel.h		\
	reader.h

# compile normally:
# make -g -O2 to optimize
all: COPT += -O2
all: CUDAOPT += -O2
all: part-pair

# compile with debug flags
debug: COPT += -DDEBUG -g
debug: CUDAOPT += -DDEBUG -g
debug: part-pair

OBJS = $(addprefix $(SRC_DIR)/, $(addsuffix .o, $(basename $(SRCC))))
OBJSCUDA = $(addprefix $(SRC_DIR)/, $(addsuffix .o, $(basename $(SRCCUDA))))

$(OBJSCUDA):$(SRC_DIR)/%.o:$(SRC_DIR)/%.cu
	$(NVCC) $(CUDAOPT) -dc $< $(LDINCS) -o $@

$(OBJS):$(SRC_DIR)/%.o:$(SRC_DIR)/%.c
	$(GCC) $(COPT) -c $< $(LDINCS) -o $@

$(SRC_DIR)/bblib.o:$(OBJSCUDA)
	$(NVCC) $(CUDAOPT) -dlink $+ $(CUDALIBS) -o $(SRC_DIR)/bblib.o

part-pair: $(OBJSCUDA) $(SRC_DIR)/bblib.o $(OBJS)
	$(GCC) $(COPT) $+ $(LDLIBS) -lcudadevrt -lstdc++ $(CUDALIBS) -o $(SIM_DIR)/part-pair

clean:
	rm -f $(SRC_DIR)/*.o $(SIM_DIR)/part-pair
