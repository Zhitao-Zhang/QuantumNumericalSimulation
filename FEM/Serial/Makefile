# =============================================================================
# 二维速度–应力 P1 三角有限元（串行 CPU）
# 参数与 Quantum/FDM/ImplicitV2/main.cpp 完全对齐。
#
# 用法：
#   make                       # 编译，用 config.h 缺省参数（200x200 x 4500 步）
#   make run                   # 编译并运行
#   make small                 # 小算例：50x50, 800 步（快速冒烟测试）
#   make clean / cleanall      # 清可执行 / 清所有 .dat
#
# 常用编译宏（同 FDM）：
#   -DNZ=200 -DNX=200 -DNT=4500  网格/时间步
#   -DDX_VAL=10 -DDZ_VAL=10      网格步长 (m)
#   -DDT_VAL=1.0e-4              时间步长 (s)
#   -DF0_HZ=30.0                 Ricker 主频 (Hz)
#   -DNBH=40                     高速夹层开关（0=均匀）
#   -DBICG_TOL=1e-9              BiCGSTAB 收敛门槛
# =============================================================================

CXX      ?= g++
CXXSTD   ?= -std=c++14
OPT      ?= -O3
WARN     ?= -Wall -Wno-unused-result
CXXFLAGS ?= $(OPT) $(WARN) $(CXXSTD)

TARGET   := fem_elastic_2d

SRCS := main.cpp \
        mesh.cpp \
        material.cpp \
        sparse_matrix.cpp \
        solver_bicgstab.cpp \
        assembler.cpp \
        field_io.cpp

OBJS := $(SRCS:.cpp=.o)

HDRS := config.h mesh.h material.h sparse_matrix.h linear_solver.h assembler.h field_io.h

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

run: $(TARGET)
	./$(TARGET)

# 快速冒烟测试：50 x 50 网格，800 步
small: CXXFLAGS += -DNZ=50 -DNX=50 -DNT=800 -DNBH=0 -DSNAP_STRIDE=20
small: clean $(TARGET)
	./$(TARGET)

.PHONY: all run small clean cleanall
clean:
	rm -f $(TARGET) $(OBJS)

cleanall: clean
	rm -f *.dat
