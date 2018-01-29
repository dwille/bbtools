## DEPENDENCIES ##
ifeq ($(USER),dwille)
  $(info User is $(USER))
  HDF5_DIR = /opt/phdf5/1.8.19
  CGNS_DIR = /opt/pcgns/3.3.0
  CUDA_DIR = /opt/cuda
  CUDA_SDK_DIR = /opt/cuda-7.5/samples
else ifeq ($(USER),dwillen3@jhu.edu)
  $(info USER is dwillen3@jhu.edu)
  HDF5_DIR = /home-1/dwillen3@jhu.edu/.deps/phdf5/1.8.18
  CGNS_DIR = /home-1/dwillen3@jhu.edu/.deps/pcgns/3.3.0
endif	

## COMPILERS ##
GCC = gcc
NVCC = nvcc

SIM_DIR = sim
SRC_DIR = src

COPT = -std=c99 -pedantic -Wall -Wextra
LDINCS = -I $(CGNS_DIR)/include # -I $(CUSP_DIR)
LDLIBS = -lm -L $(HDF5_DIR)/lib -L $(CGNS_DIR)/lib -lcgns -lhdf5

CUDAOPT = -arch=sm_30 -Xcompiler -m64

CUDAINCS = -I $(CUDA_SDK_DIR)/common/inc
CUDALIBS = -L $(CUDA_DIR)/lib64 -lcudart

SRCC =	main.c	\
	reader.c

SRCCUDA = 	cuda_fun.cu	\
		fun_kernel.cu
		

EXTRA =	Makefile	\
	main.h		\
	reader.h	\

# compile normally:
# make -g -O2 to optimize
all: COPT += -O2
all: CUDAOPT += -O2
all: single-particle-dispersion

# compile with debug flags
debug: COPT += -DDEBUG -g
debug: CUDAOPT += -DDEBUG -g
debug: single-particle-dispersion

OBJS = $(addprefix $(SRC_DIR)/, $(addsuffix .o, $(basename $(SRCC))))
OBJSCUDA = $(addprefix $(SRC_DIR)/, $(addsuffix .o, $(basename $(SRCCUDA))))

$(OBJSCUDA):$(SRC_DIR)/%.o:$(SRC_DIR)/%.cu
	$(NVCC) $(CUDAOPT) -dc $< $(CUDAINCS) $(LDINCS) -o $@

$(OBJS):$(SRC_DIR)/%.o:$(SRC_DIR)/%.c
	$(GCC) $(COPT) -c $< $(LDINCS) -o $@

$(SRC_DIR)/bblib.o:$(OBJSCUDA)
	$(NVCC) $(CUDAOPT) -dlink $+ $(CUDALIBS) -o $(SRC_DIR)/bblib.o


##  single-particle-dispersion: $(OBJS)
##  	$(GCC) $(COPT) $+ $(LDLIBS) -o $(SIM_DIR)/single-particle-dispersion
single-particle-dispersion: $(OBJSCUDA) $(SRC_DIR)/bblib.o $(OBJS)
	$(GCC) $(COPT) $+ $(LDLIBS) -lcudadevrt -lstdc++ $(CUDALIBS) -o $(SIM_DIR)/single-particle-dispersion

clean:
	rm -f $(SRC_DIR)/*.o $(SIM_DIR)/single-particle-dispersion
