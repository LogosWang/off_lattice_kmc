CXX = clang++

# 定义编译器选项
# -std=c++11：使用 C++11 标准
# -Wall：开启所有常见的警告
# -O3：开启最高级别优化
# -g：包含调试信息 (如果需要调试，可以加上；发布版本可以去掉)
CXXFLAGS = -std=c++11 -Wall -O3 -g

# 定义你的所有源文件
SRCS = main.cpp KMC_Simulator.cpp Site.cpp Event.cpp StressPoint.cpp

# 定义最终的可执行文件名称
TARGET = kmc_sim

# 默认目标：当只输入 'make' 时执行
all: $(TARGET)

# 如何从 .o 文件链接生成可执行文件
$(TARGET): $(SRCS:.cpp=.o)
	$(CXX) $(CXXFLAGS) $^ -o $@

# 如何从 .cpp 文件编译生成 .o 文件（通用规则）
# $< 代表第一个依赖文件 (.cpp 文件)
# $@ 代表目标文件 (.o 文件)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 清理目标：删除编译生成的文件
.PHONY: clean
clean:
	rm -f $(TARGET) *.o