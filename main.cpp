#include "KMC_Simulator.h" // 包含 KMC_Simulator 类的定义，这样我们才能创建和使用它
#include <iostream>        // 用于标准输出，比如打印一些程序启动信息
#include <ctime>           // 用于获取当前时间，作为随机数生成器的种子，确保每次运行结果不同
const char* KMC_SIM_CHECK = "KMC_Simulator class is visible.";

int main() {
    // --- 1. 定义模拟的基本参数 ---
    // 这些参数决定了你的模拟世界的规模和运行方式
    std::cout << KMC_SIM_CHECK << std::endl; // 你甚至可以尝试打印它
    int num_sites = 500;      // 模拟中的粒子（Site）总数量。你可以根据需要调整。
    double box_size = 5e-4;        // 模拟盒子的边长。假设是一个边长为 10.0 单位的立方体。
    double max_sim_time = 100000.0;  // 模拟将运行的最大总时间。当模拟时间达到这个值时，程序会停止。
    long long int max_sim_steps = 100000000000;    
    double jump_distance = 5e-6; 
    int CSVflag=1;
    // --- 2. 设置随机数生成器的种子 ---
    // 使用当前系统时间作为种子，可以确保每次运行程序时，随机数序列都不同，
    // 从而得到不同的模拟轨迹。如果你需要复现某个特定的模拟结果，可以使用一个固定的整数作为种子。
    unsigned int seed = static_cast<unsigned int>(time(0));

    // --- 3. 创建 KMC_Simulator 对象 ---
    // 使用我们定义的参数来实例化 KMC_Simulator。
    // 构造函数 KMC_Simulator(num_sites, box_size, seed, jump_distance) 会被调用。
    KMC_Simulator simulator(num_sites, box_size, seed, jump_distance, CSVflag);

    // --- 4. 初始化站点 ---
    // 调用 KMC_Simulator 的 initialize_sites 方法，将粒子随机放置在模拟盒子内。
    // 这会在模拟开始前设置好系统的初始状态。
    simulator.initialize_sites();

    // --- 5. 运行 KMC 模拟 ---
    // 启动 KMC 模拟的核心循环。
    // 模拟将持续进行，直到达到 `max_sim_time` 或 `max_sim_steps` 中的任意一个条件。
    simulator.run(max_sim_time, max_sim_steps);

    // --- 6. 程序结束 ---
    // 如果一切顺利，程序会返回 0，表示成功执行。
    return 0;
}