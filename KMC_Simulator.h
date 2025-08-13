#ifndef KMC_SIMULATOR_H // 防止头文件被重复包含
#define KMC_SIMULATOR_H

// --- 必需的 C++ 标准库头文件 ---
#include <vector>     // 用于 std::vector 来存储 Site 和 Event 对象
#include <string>     // 可能用于文件名或状态字符串
#include <numeric>    // 用于 std::accumulate 函数，计算总倾向
#include <iostream>   // 用于标准输入输出，比如打印模拟状态
#include <fstream>    // 用于文件操作，比如将 Site 坐标写入文件
#include <cmath>      // 用于数学函数，比如 std::log (吉莱斯皮算法中计算时间步长)
#include <iomanip>    // 用于格式化输出，比如设置浮点数精度
#include <random>     // 用于 C++11 的随机数生成器 (std::mt19937, std::uniform_real_distribution)
#include <unordered_map> // 用于 std::unordered_map，实现站点到事件索引的高效映射

// --- 你自定义的类头文件 ---
#include "Site.h"     // KMC_Simulator 需要知道 Site 类的定义来管理 Site 对象
#include "Event.h"    // KMC_Simulator 需要知道 Event 类的定义来管理 Event 对象
#include "StressPoint.h" 
class KMC_Simulator {
public:
   
    KMC_Simulator(int num_sites, double box_size, unsigned int seed, double unit_jump_distance, int CSVflag);

    
    ~KMC_Simulator();

    void initialize_sites();
    void run(double max_time, long long int max_steps);
    bool initialize_sites_from_csv(const std::string& filename);
    bool read_stress_field_from_csv(const std::string& filename);
private:
    
    static const double KB;                 // 玻尔兹曼常数 (J/K)
    static const double TEMPERATURE;        // 绝对温度 (K)
    static const double MIGRATION_BARRIER_EV; // 迁移能垒 (eV)
    static const double EV_TO_JOULE;        // eV 到焦耳的转换因子
    static const double ACTIVATION_VOLUME;  // 活化体积 (m^3) - 确保此值为正
    static const double PRE_FACTOR_V0;      // 前因子 (s^-1)


    int CSVflag;
    int num_sites;     // 粒子数量
    double box_size;   // 模拟盒子的尺寸
    double unit_jump_distance;
    double jump_distance; // 随机游走事件的跳跃步长，作为模拟器的属性

    std::mt19937 generator;         // Mersenne Twister 伪随机数引擎，提供高质量的随机数
    std::uniform_real_distribution<double> distribution; // 定义一个在 [0.0, 1.0) 范围内均匀分布的随机数分发器

    std::vector<StressPoint> stress_field_data;
    std::vector<Site> sites;       // 存储所有 Site 对象的列表，这是我们模拟的物理系统
    std::vector<Event> all_events; // 存储当前所有可能发生的具体 Event 实例的列表。
                                   // 每个 Event 对象包含其类型、倾向和关联的 Site 引用。

    // 与 all_events 列表并行，存储每个事件的倾向值。
    // 这是吉莱斯皮算法在选择事件时需要的所有倾向值的汇总。
    std::vector<double> event_propensities_for_selection; 

    // 用于高效查找的映射：将 Site 的 ID 映射到所有与该 Site 相关的 Event 在 all_events 列表中的索引。
    // 这对于局部更新事件倾向非常关键。
    std::unordered_map<int, std::vector<int>> site_to_event_indices; 
    double current_time;           // 模拟已经进行的总时间
    long long total_steps;         // 模拟已经完成的 KMC 步骤总数

    // --- KMC 核心算法的内部辅助函数 ---

    // **关键函数：计算并更新所有事件的倾向，同时构建 all_events 列表**
    // 每次模拟状态改变（如一个事件发生）后，都需要调用此函数以更新事件池。
    int find_closest_stress_point_id(double px, double py, double pz) const;
    void calculate_all_propensities_and_events(); 
    // **关键函数：根据吉莱斯皮算法选择下一个要执行的事件**
    // 返回被选中的 Event 在 all_events 列表中的索引。
    // 如果没有事件可发生，则返回 -1。
    int select_event_index();

    void execute_event(int event_index); 

    // 应用周期性边界条件 (PBC)：确保粒子在模拟盒子内。
    void apply_pbc(Site& s);

    // **关键优化函数：在事件发生后，只更新受影响站点（及其邻居）的事件倾向**
    // affected_site_id: 刚刚发生事件的 Site 的 ID。
    double calculate_site_random_walk_propensity(const Site& s) const;
    void update_affected_events(int affected_site_id);

    double get_uniform_random();      // 生成一个 [0.0, 1.0) 范围内的均匀随机实数
    int get_uniform_int(int min, int max); // 生成一个在 [min, max] 范围内的均匀随机整数

    // --- 模拟状态和输出函数 ---
    void print_status(); // 打印当前模拟时间、步数等状态信息
    void dump_sites(long long int step); // 将所有 Site 的当前坐标输出到文件，用于后续分析或可视化
};

#endif // KMC_SIMULATOR_Hs