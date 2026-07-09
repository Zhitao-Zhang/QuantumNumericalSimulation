# =============================================================================
# 二维位移二阶弹性波 P1 三角有限元（串行 CPU） — 方案 B / SESSION 7
# 参数与 Quantum/FDM/ImplicitV2/main.cpp 对齐。
#
# 格式：M ü + C u̇ + K u = f，Newmark 平均加速度隐式（无条件稳定）。
# 每步解 A a = b（A = M + γ·dt·C + β·dt²·K），后端 Jacobi-BiCGSTAB。
#
# 用法：
#   make                       # 200x200 x 4100 步，默认参数（层状介质 + 吸收层）
#   make run                   # 编译并运行
#   make homo                  # 均匀介质 (NBH=0)，检验波前是否为干净圆形
#   make small                 # 50x50 x 800 步，快速冒烟测试
#   make clean / cleanall      # 清可执行 / 清所有 .dat
#
# 常用编译宏：
#   -DNZ=200 -DNX=200 -DNT=4100  网格/时间步
#   -DDX_VAL=10 -DDZ_VAL=10      网格步长 (m)
#   -DDT_VAL=1.0e-4              时间步长 (s)
#   -DF0_HZ=30.0                 Ricker 主频 (Hz)
#   -DNBH=40                     高速夹层开关（0=均匀）
#   -DNABSDEF=30                 吸收层厚度 cells（0=关闭）
#   -DNEWMARK_BETA / -DNEWMARK_GAMMA   Newmark 参数
#   -DSRC_SCALE=...             声源矩幅度（默认 dx*dz）
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

# 均匀介质：最干净地检验波前是否为圆形（无层界面干扰）
homo: CXXFLAGS += -DNBH=0
homo: clean $(TARGET)
	./$(TARGET)

# 快速冒烟测试：50 x 50 网格，800 步
small: CXXFLAGS += -DNZ=50 -DNX=50 -DNT=800 -DNBH=0 -DSNAP_STRIDE=20 -DNABSDEF=10
small: clean $(TARGET)
	./$(TARGET)

# 均匀介质小算例（快速看波前形状）
small_homo: CXXFLAGS += -DNZ=120 -DNX=120 -DNT=700 -DNBH=0 -DSNAP_STRIDE=20 -DNABSDEF=15
small_homo: clean $(TARGET)
	./$(TARGET)

.PHONY: all run homo small small_homo clean cleanall
clean:
	rm -f $(TARGET) $(OBJS)

cleanall: clean
	rm -f *.dat
