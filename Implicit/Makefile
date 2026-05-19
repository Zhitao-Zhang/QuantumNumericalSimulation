# ==== Toolchains ====
CXX        := mpic++
NVCC       := nvcc -ccbin=$(CXX)

# ==== Paths (按你当前 Spack 环境，可按需覆盖 CUDA_HOME) ====
CUDA_HOME  := /home/spack/spack/opt/spack/linux-debian12-cascadelake/gcc-12.2.0/cuda-12.8.0-ais4c4h3sfikeullxq6djwouyp72pfvd
MPI_HOME   := /home/spack/spack/opt/spack/linux-debian12-cascadelake/gcc-12.2.0/openmpi-5.0.6-r2esph5hsmnvadmhvjvngfiqx2xqlhit

INC_DIRS   := -I$(CUDA_HOME)/include -I$(MPI_HOME)/include
LIB_DIRS   := -L$(CUDA_HOME)/lib64
RPATH_DIRS := $(CUDA_HOME)/lib64

# ==== Architectures ====
# A100 固定 80。把 5090D �? SM 改成你的实际值（示例�? 89）：
SM_A100   := 80
SM_5090D  := 89

GENCODE_FLAGS := \
  -gencode arch=compute_$(SM_A100),code=sm_$(SM_A100) \
  -gencode arch=compute_$(SM_A100),code=compute_$(SM_A100) \
  -gencode arch=compute_$(SM_5090D),code=sm_$(SM_5090D) \
  -gencode arch=compute_$(SM_5090D),code=compute_$(SM_5090D)

# ==== Flags ====
CXXFLAGS   := -O3 -w $(INC_DIRS)
NVCCFLAGS  := -O3 -Xcompiler -fPIC -lineinfo $(INC_DIRS) $(GENCODE_FLAGS)
LDFLAGS    := $(LIB_DIRS) -Xlinker -rpath -Xlinker $(RPATH_DIRS)
LDLIBS     := -lcudart

# ==== Targets ====
EXECUTABLE := my_program
SRC_CU     := main.cu
OBJ        := $(SRC_CU:.cu=.o)

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJ)
	$(NVCC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

# 通用 .cu �? .o 规则
%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJ) $(EXECUTABLE)