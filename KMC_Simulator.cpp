#include "KMC_Simulator.h" // 包含 KMC_Simulator 的头文件
#include <numeric>         // 用于 std::accumulate
#include <algorithm>       // 可能用于未来更复杂的事件更新逻辑 (如 std::remove_if)
#include <iostream>        // 用于调试输出和状态报告
#include <sstream>
#include <limits>

// --- 静态常量成员的定义和初始化 ---
// 这些值应该放在 KMC_Simulator 类的成员函数之外
// **** 请根据 Si 在奥氏体不锈钢中的实际物理参数调整这些值 ****
const double KMC_Simulator::KB = 1.380649e-23; // J/K (玻尔兹曼常数)
const double KMC_Simulator::TEMPERATURE = 800.0 + 273.15; // K (例如 800 摄氏度，转换为开尔文)
const double KMC_Simulator::MIGRATION_BARRIER_EV = 1.0; // eV (示例迁移能垒，请查阅文献)
const double KMC_Simulator::EV_TO_JOULE = 1.60218e-19; // J/eV (eV 到焦耳的转换因子)
const double KMC_Simulator::ACTIVATION_VOLUME = 1.0e-29; // m^3 (示例活化体积，请查阅文献，通常为正值)
const double KMC_Simulator::PRE_FACTOR_V0 = 1.0e13; // s^-1 (示例前因子，请查阅文献)


KMC_Simulator::KMC_Simulator(int num_sites_arg, double box_size_arg, unsigned int seed, double unit_jump_distance_arg, int CSVflag) : 
    num_sites(num_sites_arg),             // 粒子数量
    box_size(box_size_arg),               // 盒子大小
    unit_jump_distance(unit_jump_distance_arg),
    CSVflag(CSVflag),     // 随机游走步长
    generator(seed),                      // 初始化随机数引擎，使用提供的种子
    distribution(0.0, 1.0),               // 初始化均匀分布器，范围 [0.0, 1.0)
    current_time(0.0),                    // 模拟开始时时间为 0
    total_steps(0)                        // 模拟开始时步数为 0
{
    // 预留 `sites` 向量的空间，避免在添加 Site 时不必要的内存重新分配
    sites.reserve(num_sites); 
    
    // all_events 和 event_propensities_for_selection 的大小会在 calculate_all_propensities_and_events 中动态调整
    std::cout << "KMC Simulator initialized with " << num_sites << " sites, box size " << box_size 
              << ", and jump distance " << jump_distance << "." << std::endl;
    read_stress_field_from_csv("stress.csv");
}

KMC_Simulator::~KMC_Simulator() {
    std::cout << "KMC Simulator shut down." << std::endl;
}
double KMC_Simulator::get_uniform_random() {
    return distribution(generator); 
}
int KMC_Simulator::get_uniform_int(int min, int max) {
    // 每次调用时创建临时的分发器，确保每次都能正确处理不同的范围
    std::uniform_int_distribution<int> int_distribution(min, max);
    return int_distribution(generator);
}
void KMC_Simulator::initialize_sites() {
    if (CSVflag==0){
    for (int i = 0; i < num_sites; ++i) {
        // 为每个 Site 生成随机的 X, Y, Z 坐标
        double x = get_uniform_random() * box_size; 
        double y = get_uniform_random() * box_size; 
        double z = get_uniform_random() * box_size; 
        
        // 使用 emplace_back 高效地在向量末尾构造 Site 对象
        sites.emplace_back(i, x, y, z, 0); // Site 构造函数: Site(id, x, y, z)
    }
    std::cout << "Sites initialized randomly within the box." << std::endl;
    print_status(); // 打印当前模拟状态
}
else {
    initialize_sites_from_csv("sites.csv");
}
}


bool KMC_Simulator::initialize_sites_from_csv(const std::string& filename) {
    std::ifstream infile(filename); // 打开文件
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open CSV file: " << filename << std::endl;
        return false; // 文件打开失败
    }

    sites.clear(); // 清空当前站点列表，准备从 CSV 读取

    std::string line;
    // 读取并跳过 CSV 文件的第一行（通常是列头）
    if (std::getline(infile, line)) { 
        std::cout << "Skipping CSV header: " << line << std::endl;
    } else {
        std::cerr << "Error: CSV file is empty or could not read header." << std::endl;
        return false;
    }

    int id;
    double x, y, z;
    int type;
    char comma; // 用于读取逗号分隔符

    // 逐行读取数据
    while (std::getline(infile, line)) {
        std::stringstream ss(line); // 将行字符串放入字符串流

        // 从字符串流中解析数据，用逗号分隔
        // 确保你的 CSV 文件格式是 id,x,y,z,type
        if (ss >> id >> comma >> x >> comma >> y >> comma >> z >> comma >> type) {
            sites.emplace_back(id, x, y, z, type);
        } else {
            std::cerr << "Warning: Skipping malformed line in CSV: " << line << std::endl;
        }
    }

    infile.close(); // 关闭文件

    // 更新 num_sites 以匹配从 CSV 文件中读取的站点数量
    num_sites = sites.size(); 
    if (num_sites == 0) {
        std::cerr << "Error: No sites were read from the CSV file." << std::endl;
        return false;
    }

    std::cout << "Sites initialized from CSV file: " << filename << ". Total sites: " << num_sites << std::endl;
    print_status();
    return true; // 成功读取
}

bool KMC_Simulator::read_stress_field_from_csv(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "错误：无法打开应力场文件 " << filename << std::endl;
        return false;
    }

    std::string line;
    // 假设第一行是标题行，跳过它
    if (std::getline(file, line)) {
        std::cout << "skip first line: " << line << std::endl;
    } else {
        std::cerr << "错误：应力场文件为空或无法读取文件头。" << std::endl;
        return false;
    }

    int current_id = 0; // 为 StressPoint 分配一个连续的 ID
    while (std::getline(file, line)) {
        if (line.empty()) continue; // 跳过空行

        std::stringstream ss(line);
        std::string segment;
        std::vector<double> values;

        // 使用逗号作为分隔符解析每一行
        while(std::getline(ss, segment, ',')) {
            try {
                values.push_back(std::stod(segment));
            } catch (const std::invalid_argument& e) {
                std::cerr << "错误：解析应力场文件时无效的数字转换：" << segment << " in line: " << line << std::endl;
                return false;
            } catch (const std::out_of_range& e) {
                std::cerr << "错误：解析应力场文件时数字超出范围：" << segment << " in line: " << line << std::endl;
                return false;
            }
        }

        // 检查是否有足够的列（X, Y, Z, Trace）
        if (values.size() >= 4) { // 至少需要4列
            // 如果 CSV 中没有 ID，我们自己生成一个
            stress_field_data.emplace_back(current_id++, values[0], values[1], values[2], values[3]);
        } else {
            std::cerr << "警告：应力场文件行格式不正确，跳过此行： " << line << std::endl;
        }
    }

    std::cout << "successfully read " << stress_field_data.size() << " stress point" << std::endl;
    return true;
}

int KMC_Simulator::find_closest_stress_point_id(double px, double py, double pz) const {
    if (stress_field_data.empty()) {
        std::cerr << "错误: 应力场数据为空，无法查找最近点。" << std::endl;
        return -1; // 或者可以考虑抛出异常
    }

    double min_dist = std::numeric_limits<double>::max(); // 存储最小距离
    int closest_id = -1; // 存储最近点的 ID

    for (const auto& sp : stress_field_data) {
        // 使用 StressPoint 类的 get_distance 方法来计算距离
        double current_dist = sp.get_distance(px, py, pz);

        if (current_dist < min_dist) {
            min_dist = current_dist;
            closest_id = sp.id;
        }
    }
    return closest_id;
}


void KMC_Simulator::calculate_all_propensities_and_events() {
    all_events.clear();                 
    event_propensities_for_selection.clear(); 
    site_to_event_indices.clear();      

    int current_event_index = 0; // 用于跟踪当前事件在 all_events 向量中的索引
    for (int i = 0; i < num_sites; ++i) {
        // **计算事件倾向：**
        // 对于随机游走事件，我们假设其倾向是一个固定值，例如 1.0。
        // 如果未来引入依赖于 Site 状态（如是否有空位）的倾向，这部分逻辑会更复杂。
        double calculated_propensity = calculate_site_random_walk_propensity(sites[i]); 

        // **创建 Event 对象：**
        // 使用 Event 构造函数 Event(EventType type, double prop, const Site& s)
        // 注意：这里仍然传入了 `sites[i]` 的引用。
        // 我再次强调，虽然你的 `Event` 类定义允许这样做，但它有潜在的风险（向量重新分配导致引用失效）。
        // 最安全的做法是 `Event` 存储 `sites[i].id`，然后在这里传入 ID。
        Event rw_event(1, calculated_propensity, sites[i]); 
        
        // 将创建的 Event 对象添加到总事件列表 `all_events`
        all_events.push_back(rw_event); 
        // 将该事件的倾向值添加到用于选择的倾向列表 `event_propensities_for_selection`
        event_propensities_for_selection.push_back(calculated_propensity); 
        
        // **更新站点-事件索引映射：**
        // 将当前 Site 的 ID 映射到 `all_events` 中它所关联的事件的索引。
        // 一个随机游走事件只与一个站点关联。
        site_to_event_indices[sites[i].id].push_back(current_event_index);
        
        current_event_index++; 
    }
}

double KMC_Simulator::calculate_site_random_walk_propensity(const Site& s) const {
    // 1. 找到距离当前 Site 最近的 StressPoint 的 ID (即 vector 索引)
    int closest_sp_id = find_closest_stress_point_id(s.x, s.y, s.z);

    double site_stress_trace = 0.0; // 初始化应力迹

    // 2. 根据 ID (索引) 从 stress_field_data (vector) 中获取对应的 StressPoint 的 trace 值
    //    这里进行了边界检查，确保 closest_sp_id 是一个有效的索引
    if (closest_sp_id >= 0 && closest_sp_id < static_cast<int>(stress_field_data.size())) {
        site_stress_trace = stress_field_data[closest_sp_id].trace;
    } else {
        // 如果 find_closest_stress_point_id 返回 -1 或 ID 超出 vector 范围，
        // 则打印警告并使用默认应力迹 0.0。
        std::cerr << "警告：未找到 Site " << s.id << " 对应的应力场数据点 (ID: " << closest_sp_id << ")，使用默认应力迹 0.0。" << std::endl;
        site_stress_trace = 0.0;
    }

    // 3. 使用你最终确认的公式计算倾向：
    //    v = v0 * exp(-(DeltaE + tr(sigma) * Omega) / (kB * T))
    //    其中 DeltaE 已经从 eV 转换为焦耳 (DELTA_E_JOULE)
    double DELTA_E_JOULE = MIGRATION_BARRIER_EV * EV_TO_JOULE;

    // 计算有效能垒： DeltaE + tr(sigma) * Omega
    double effective_barrier = DELTA_E_JOULE + (site_stress_trace * ACTIVATION_VOLUME);

    // 计算指数项的除数： kB * T
    double kb_T = KB * TEMPERATURE;

    // 计算完整的指数项： effective_barrier / (kB * T)
    double exponent_term = effective_barrier / kb_T;

    // 最终计算倾向： v0 * exp(-exponent_term)
    double propensity = PRE_FACTOR_V0 * std::exp(-exponent_term);
    propensity *= 1e-3; 
    // ms-1 unit

    // 4. 确保计算出的倾向为正数，防止数学错误或不合理结果
    if (propensity <= 0) {
        propensity = 1e-30; // 设置一个非常小的正值，避免零或负倾向
        std::cerr << "警告：为 Site " << s.id << " 计算出的倾向为非正值 (" << propensity << ")，已强制设置为 1e-30。最近应力迹: " << site_stress_trace << std::endl;
    }

    return propensity;
}


int KMC_Simulator::select_event_index() {
    // 计算所有事件的总倾向 H (吉莱斯皮算法中的 A_0)
    double total_propensity = std::accumulate(event_propensities_for_selection.begin(), event_propensities_for_selection.end(), 0.0);

    // 如果总倾向为零或负，说明没有事件可以发生
    if (total_propensity <= 0.0) {
        return -1;
    }

    // 生成一个随机数 r1，用于选择事件。r1 * H 将在 [0, H) 范围内
    double r1_times_H = get_uniform_random() * total_propensity; 
    
    double sum_partial_propensity = 0.0;
    int chosen_event_index = -1;

    // 遍历事件倾向列表，累加倾向值，直到累加和大于 r1_times_H
    // 第一个满足条件的事件就是被选中的事件
    for (int i = 0; i < all_events.size(); ++i) { 
        sum_partial_propensity += event_propensities_for_selection[i];
        if (sum_partial_propensity > r1_times_H) {
            chosen_event_index = i; // 找到了被选中的事件的索引
            break;
        }
    }
    return chosen_event_index;
}

void KMC_Simulator::execute_event(int event_index) {
    // 获取被选中的 Event 对象的常量引用
    const Event& chosen_event = all_events[event_index];
    

    Site& target_site = const_cast<Site&>(chosen_event.site); 

    // 根据 Event 对象的类型 (etype) 来执行不同的操作
    // KMC_Simulator 现在承担了所有事件类型的执行逻辑。
    switch (chosen_event.etype) {
        case 1: {
            // 实现随机游走事件的逻辑：
            std::uniform_int_distribution<int> axis_dist(0, 2);   // 选择轴 (X, Y, Z)
            std::uniform_int_distribution<int> sign_dist(0, 1);   // 选择方向 (+/-)

            int axis = axis_dist(generator);                  // 随机选择移动的轴
            double direction_sign = (sign_dist(generator) == 0 ? -1.0 : 1.0); // 随机选择移动的方向
            
            double displacement = jump_distance * direction_sign;

            // 调用 Site 对象的 `move` 方法来更新其位置
            target_site.move(axis, displacement);
            apply_pbc(target_site); 
            double updated_prop=calculate_site_random_walk_propensity(target_site);
            all_events[event_index].propensity=updated_prop;
            event_propensities_for_selection[event_index]=updated_prop;
            break;
        }
        default: {
            // 如果遇到未知的事件类型，则输出错误信息
            std::cerr << "Error: Unknown event type (" << chosen_event.etype 
                      << ") encountered in KMC_Simulator::execute_event!" << std::endl;
            break;
        }
    }
    
    // 事件执行后，对受影响的 Site 应用周期性边界条件，确保它仍在模拟盒子内
    // apply_pbc(target_site); 
}

void KMC_Simulator::apply_pbc(Site& s) {
    // `fmod` 计算浮点数余数，确保在 [0, box_size) 范围内
    s.x = fmod(s.x, box_size); 
    if (s.x < 0) s.x += box_size; // 如果结果为负，加上 box_size 使其为正

    s.y = fmod(s.y, box_size); 
    if (s.y < 0) s.y += box_size;
    
    s.z = fmod(s.z, box_size); 
    if (s.z < 0) s.z += box_size;
}

void KMC_Simulator::update_affected_events(int affected_site_id) {
    // // 检查 `affected_site_id` 是否在映射中存在
    // if (site_to_event_indices.count(affected_site_id)) { 
    //     // 获取所有与 `affected_site_id` 关联的事件的索引列表
    //     const auto& related_event_indices = site_to_event_indices.at(affected_site_id);

    //     for (int event_idx : related_event_indices) {
    //         // 获取待更新事件的引用
    //         Event& event_to_update = all_events[event_idx];
            
    //         double new_propensity = 0.0;
    //         // 根据事件类型重新计算其倾向
    //         switch (event_to_update.etype) {
    //             case RANDOM_WALK:
    //                 // 对于随机游走，我们假设倾向始终为 1.0
    //                 new_propensity = 1.0; 
    //                 // 如果未来倾向依赖于 Site 的状态（如是否空），这里需要访问 Site 对象
    //                 // 例如：Site& affected_site = const_cast<Site&>(event_to_update.site);
    //                 // if (affected_site.is_empty()) new_propensity = 1.0; else new_propensity = 0.0;
    //                 break;
    //             // --- 未来扩展点：其他事件类型的倾向计算 ---
    //             /*
    //             case ADSORPTION:
    //                 // new_propensity = calculate_adsorption_propensity(const_cast<Site&>(event_to_update.site));
    //                 break;
    //             */
    //             default:
    //                 std::cerr << "Error: Unhandled event type for propensity update." << std::endl;
    //                 break;
    //         }
    //         // 更新 Event 对象内部存储的倾向值
    //         event_to_update.propensity = new_propensity;
    //         // **最重要的是，更新用于选择事件的倾向列表！**
    //         event_propensities_for_selection[event_idx] = new_propensity;
    //     }
}

void KMC_Simulator::print_status() {
    std::cout << "Time: " << std::fixed << std::setprecision(6) << current_time 
              << " | Steps: " << total_steps << std::endl;
}

// --- 输出 Site 坐标到文件 ---
// 将所有 Site 的当前三维坐标以 .xyz 格式写入文件，方便外部工具（如 VMD）进行可视化。
void KMC_Simulator::dump_sites(int step) {
    std::string filename = "sites_" + std::to_string(step) + ".xyz"; // 文件名包含 KMC 步数
    std::ofstream outfile(filename); // 创建并打开文件

    if (!outfile.is_open()) { // 检查文件是否成功打开
        std::cerr << "Error: Could not open file " << filename << " for dumping sites." << std::endl;
        return;
    }

    // .xyz 文件格式要求：
    outfile << num_sites << std::endl; // 第一行：原子数量
    outfile << "KMC Step: " << step << ", Time: " << std::fixed << std::setprecision(6) << current_time << std::endl; // 第二行：注释

    // 写入每个 Site 的数据
    for (const auto& s : sites) {
        outfile << "Si " << std::fixed << std::setprecision(6) 
                << s.x << " " << s.y << " " << s.z << std::endl;
        // std::cout<<"writing"<<std::endl;
    }
    outfile.close(); // 关闭文件
}

void KMC_Simulator::run(double max_time, int max_steps) {
    std::cout << "\nStarting KMC simulation..." << std::endl;
    std::cout << "Maximum Simulation Time: " << max_time 
              << ", Maximum KMC Steps: " << max_steps << std::endl;

    // **首次构建事件列表和计算所有倾向**
    calculate_all_propensities_and_events();
    dump_sites(0); // 输出初始构型

    // KMC 模拟主循环
    while (current_time < max_time && total_steps < max_steps) {
        // 获取当前所有事件的总倾向 H
        double total_propensity = std::accumulate(event_propensities_for_selection.begin(), event_propensities_for_selection.end(), 0.0);

        // 如果总倾向为零，说明没有事件可以发生，模拟结束
        if (total_propensity <= 0.0) {
            std::cout << "Total propensity is zero. No more events possible. Exiting simulation." << std::endl;
            break;
        }

        // **吉莱斯皮算法：时间推进**
        double r2 = get_uniform_random(); // 生成第二个随机数 r2
        double dt = -std::log(r2) / total_propensity; // 根据 r2 和总倾向计算时间步长
        
        current_time += dt; // 更新模拟时间

        // 检查是否已达到最大模拟时间，如果达到则提前退出循环
        if (current_time >= max_time) {
            current_time = max_time; // 确保最终时间不超过设定值
            std::cout << "Reached maximum simulation time. Exiting KMC loop." << std::endl;
            break;
        }

        // **吉莱斯皮算法：事件选择**
        int chosen_event_idx = select_event_index(); // 选择下一个要执行的事件的索引

        if (chosen_event_idx == -1) {
             std::cerr << "Error: No event selected despite positive total propensity. Exiting." << std::endl;
             break;
        }

        // **执行选中的事件**
        jump_distance=unit_jump_distance*std::sqrt(dt*all_events[chosen_event_idx].propensity);
        execute_event(chosen_event_idx); 
        total_steps++; // 增加已完成的 KMC 步数

        // **事件发生后，更新事件倾向**
        // 对于你当前只有随机游走（倾向固定为 1.0）的模拟，
        // `update_affected_events` 函数可能显得多余，因为倾向不会改变。
        // 但在更复杂的模拟中，这里至关重要，它决定了 KMC 算法的正确性和性能。
        //
        // 鉴于你目前的 `Event` 类依然包含 `const Site& site` 引用，
        // 每次 KMC 步都**重新计算所有事件并重建 `all_events` 列表**
        // (`calculate_all_propensities_and_events()` 的作用) 
        // 是一个安全但效率较低的方式来确保所有引用都指向有效内存。
        //
        // 如果你已改为 `int site_id`，那么就可以调用更高效的 `update_affected_events()`：
        // update_affected_events(all_events[chosen_event_idx].site_id);
        // 并在此基础上实现局部更新。
        
       // calculate_all_propensities_and_events(); // 目前继续全量重建事件列表

        // 每隔一定步数输出模拟状态和粒子构型
        if (total_steps % 10000000 == 0) { // 例如，每 1000 步输出一次
            print_status();
            dump_sites(total_steps);
        }
    }
    // 模拟结束时，输出最终状态和构型
    std::cout << "\nSimulation finished. Final status:" << std::endl;
    print_status();
    dump_sites(total_steps); 
}