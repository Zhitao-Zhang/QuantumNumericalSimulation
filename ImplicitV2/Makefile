# 串行 CPU 版本：仅依赖 g++，无 CUDA / MPI。
# 可用宏覆盖参数：
#   -DNZ=200 -DNX=200 -DNT=2500  网格尺寸 / 时间步数
#   -DNPMLDEF=20                  PML 厚度（cells）
#   -DDX_VAL=1.0 -DDZ_VAL=1.0     x, z 方向网格步长 (m)
#   -DDT_VAL=1.0e-4               时间步长 (s)
#   -DF0_HZ=75.0                  Ricker 子波主频 (Hz)

CXX      ?= g++
CXXSTD   ?= -std=c++14
OPT      ?= -O2
WARN     ?= -Wall -Wno-unused-result
CXXFLAGS ?= $(OPT) $(WARN) $(CXXSTD)

TARGET   := cn_elastic_2d
SRC      := main.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

run: $(TARGET)
	./$(TARGET)

.PHONY: clean cleanall
clean:
	rm -f $(TARGET)

cleanall: clean
	rm -f *.dat
