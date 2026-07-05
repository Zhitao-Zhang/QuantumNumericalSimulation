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

# 快速冒烟测试：50 x 50 网格，800 步（含 PML）
small: CXXFLAGS += -DNZ=50 -DNX=50 -DNT=800 -DNBH=0 -DSNAP_STRIDE=20 -DNPMLDEF=10
small: clean $(TARGET)
	./$(TARGET)

# 无 PML 冒烟测试：验证 NPMLDEF=0 时退化为旧行为
small_nopml: CXXFLAGS += -DNZ=50 -DNX=50 -DNT=800 -DNBH=0 -DSNAP_STRIDE=20 -DNPMLDEF=0
small_nopml: clean $(TARGET)
	./$(TARGET)

# ================= Q1 四边形网格测试 =================
# 用 Q1 双线性四边形单元代替 P1 三角形（编译期 -DUSE_Q1_MESH）
# 目的：诊断 P1-P1 velocity-stress 是否是散射波包/环的根源
# 每 cell 一个 4-node 单元；每节点连接 3×3 stencil（8 邻居），
# 更接近 FDM 的中心差分行为。
#
# 用法：
#   make small_q1        # 50x50, 800 步 (Q1 + PML)  快速烟测
#   make q1              # 全尺寸 200x200 (Q1)
#   make q1_homo         # 200x200 (Q1) 均匀介质 (NBH=0)
q1: CXXFLAGS += -DUSE_Q1_MESH
q1: clean $(TARGET)
	./$(TARGET)

q1_homo: CXXFLAGS += -DUSE_Q1_MESH -DNBH=0
q1_homo: clean $(TARGET)
	./$(TARGET)

small_q1: CXXFLAGS += -DUSE_Q1_MESH -DNZ=50 -DNX=50 -DNT=800 -DNBH=0 -DSNAP_STRIDE=20 -DNPMLDEF=10
small_q1: clean $(TARGET)
	./$(TARGET)

# 均匀介质对比目标：用 P1 三角、均匀介质，能干净地检验散射是否来自离散化本身
homo: CXXFLAGS += -DNBH=0
homo: clean $(TARGET)
	./$(TARGET)

.PHONY: all run small small_nopml small_q1 q1 q1_homo homo clean cleanall
clean:
	rm -f $(TARGET) $(OBJS)

cleanall: clean
	rm -f *.dat
