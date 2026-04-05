#include "KMC_Simulator.h" // 包含 KMC_Simulator 的头文件
#include <numeric>         // 用于 std::accumulate
#include <algorithm>       // 可能用于未来更复杂的事件更新逻辑 (如 std::remove_if)
#include <iostream>        // 用于调试输出和状态报告
#include <sstream>
#include <limits>
#include <filesystem>
namespace fs = std::filesystem;

// --- 静态常量成员的定义和初始化 ---
// 这些值应该放在 KMC_Simulator 类的成员函数之外
// **** 请根据 Si 在奥氏体不锈钢中的实际物理参数调整这些值 ****
const double KMC_Simulator::KB = 1.380649e-23; // J/K (玻尔兹曼常数)
const double KMC_Simulator::TEMPERATURE = 400.0 + 273.15; // K (例如 800 摄氏度，转换为开尔文)
const double KMC_Simulator::MIGRATION_BARRIER_EV = 1.0; // eV (示例迁移能垒，请查阅文献)
const double KMC_Simulator::EV_TO_JOULE = 1.60218e-19; // J/eV (eV 到焦耳的转换因子)
const double KMC_Simulator::ACTIVATION_VOLUME = 1.0e-28; // m^3 (示例活化体积，请查阅文献，通常为正值)
const double KMC_Simulator::PRE_FACTOR_V0 = 1.0e13; // s^-1 (示例前因子，请查阅文献)
const double KMC_Simulator::OXIDATION_THRESH = 2.5e-7; 
const double KMC_Simulator::CR_OXIDATION_GIBBS = -1.1713e-18/2.0; 
const double KMC_Simulator::SI_OXIDATION_GIBBS = -1.4215e-18/2.0; 
const double KMC_Simulator::NI_OXIDATION_GIBBS = -7.0309e-19/2.0;
const double KMC_Simulator::NI_OXIDE_THICKNESS = 1.62e-3;
const double KMC_Simulator::SI_OXIDE_THICKNESS = 3.29e-3;
const double KMC_Simulator::CR_OXIDE_THICKNESS = 2.11e-3;    
const double KMC_Simulator::ALPHA = 0.02;
const double KMC_Simulator::OXIDATION_BARRIER = 3e-20; 
const int KMC_Simulator::OXYGEN_DIFF=3;
const int KMC_Simulator::SILI_DIFF=1;
const int KMC_Simulator::CROM_DIFF=4;
const int KMC_Simulator::NIK_DIFF=5;
const int KMC_Simulator::OXYGEN_ABS=2;
const int KMC_Simulator::OXYG=2;
const int KMC_Simulator::CROM=1;
const int KMC_Simulator::SILIC=0;
const int KMC_Simulator::NIK=3;
const int KMC_Simulator::RANDDIFF=0;
const int KMC_Simulator::TRACKDIFF=1;


KMC_Simulator::KMC_Simulator(int num_sites_arg, double box_size_arg, unsigned int seed, double unit_jump_distance_arg, int CSVflag, double GB_A, double GB_B, double GB_C, double GB_D, int Diff_mod) : 
    num_sites(num_sites_arg),             // 粒子数量
    box_size(box_size_arg),               // 盒子大小
    unit_jump_distance(unit_jump_distance_arg),
    DOSE(0.0),
    GB_normal_stress(0.0),
    CSVflag(CSVflag),     
    GB_A(GB_A),
    GB_B(GB_B),
    GB_C(GB_C),
    GB_D(GB_D),
    Diff_mod(Diff_mod),
    
    generator(seed),                      // 初始化随机数引擎，使用提供的种子
    distribution(0.0, 1.0),               // 初始化均匀分布器，范围 [0.0, 1.0)
    current_time(0.0),                    // 模拟开始时时间为 0
    total_steps(0),
    rate_scale(0.0001),
    num_si_oxide(0),
    num_cr_oxide(0),
    num_ni_oxide(0),
    num_total_oxide(0),                        // 模拟开始时步数为 0
    grain_boundary(GB_A, GB_B, GB_C, GB_D, box_size_arg),
    total_rate(0.0),
    output_dir("00"),
    cells(),
    nx(100),
    ny(100),
    nz(100),
    cell_size_x(0.001),
    cell_size_y(0.001),
    cell_size_z(0.001)
{
    // 预留 `sites` 向量的空间，避免在添加 Site 时不必要的内存重新分配
    sites.reserve(num_sites); 
    
    // all_events 和 event_propensities_for_selection 的大小会在 calculate_all_propensities_and_events 中动态调整
    std::cout << "KMC Simulator initialized with " << num_sites << " sites, box size " << box_size 
              << ", and jump distance " << jump_distance << "." << std::endl;
    read_stress_field_from_csv("stress.csv");
    // output_dir=std::to_string(static_cast<int>(10.0*DOSE));
}

KMC_Simulator::~KMC_Simulator() {
    std::cout << "KMC Simulator shut down." << std::endl;
}
double KMC_Simulator::get_uniform_random() {
    return distribution(generator); 
}

void KMC_Simulator::read_input_file(const std::string& filename)
{
    if (filename.empty()) return;  // 没指定 -i 就保持默认

    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "Warning: cannot open input file: " << filename
                  << " (using default parameters)\n";
        return;
    }

    std::string line;
    while (std::getline(fin, line)) {
        // 去注释
        auto pos = line.find('#');
        if (pos != std::string::npos) line = line.substr(0, pos);

        std::istringstream iss(line);
        std::string key;
        if (!(iss >> key)) continue;

        if (key == "dose") {
            iss >> DOSE;                 // 覆盖默认 DOSE
        } else if (key == "output_dir") {
            iss >> output_dir;           // 覆盖默认 output_dir
        } 
        else if (key == "GB_normal_stress") {
            iss >> GB_normal_stress;     // 覆盖默认 GB_normal_stress
        }
        else {
            // 最小改动：先忽略未知 key，避免影响你现有功能
            // std::cerr << "Warning: unknown key: " << key << "\n";
        }
    }

    // 目录确保存在（你已经实现过目录创建）
    if (!output_dir.empty()) {
        fs::create_directories(output_dir);
    }
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
    initialize_sites_from_csv("sites3type_more.csv");
}
}

void KMC_Simulator::init_cell_list(int count_x, int count_y, int count_z) {
    this->nx = count_x;
    this->ny = count_y;
    this->nz = count_z;

    // 1. 确定边界 (假设从 0 开始)
    double min_x =0,  min_y = 0,  min_z = 0.0; 
    double max_x = box_size, max_y = box_size, max_z = box_size;
    for (const auto& s : sites) {
        max_x = std::max(max_x, s.x);
        max_y = std::max(max_y, s.y);
        max_z = std::max(max_z, s.z);
    }

    this->cell_size_x = (max_x - min_x) / nx;
    this->cell_size_y = (max_y - min_y) / ny;
    this->cell_size_z = (max_z - min_z) / nz;
    
    // 如果三个方向格子大小一样，可以取平均或最大值，但通常建议分别存储
    // 为了简单，我们假设 cell_size 在各个方向是一样的，或者你定义三个变量

    // 3. 分配内存
    cells.assign(nx * ny * nz, std::vector<int>());

    // 4. 填入数据
    for (int i = 0; i < (int)sites.size(); ++i) {
        int idx = get_cell_index(sites[i].x, sites[i].y, sites[i].z);
        cells[idx].push_back(i);
    }
    std::cout<<"cells number: "<<cells.size()<<std::endl;
}

void KMC_Simulator::update_site_cell(int site_id, double old_x, double old_y, double old_z) {
    int old_idx = get_cell_index(old_x, old_y, old_z);
    int new_idx = get_cell_index(sites[site_id].x, sites[site_id].y, sites[site_id].z);

    // 只有在发生跨格子移动时才进行 vector 操作
    if (old_idx != new_idx) {
        auto& old_vec = cells[old_idx];
        // 这一步是从旧格子的 vector 中把当前原子的 ID 删掉
        old_vec.erase(std::remove(old_vec.begin(), old_vec.end(), site_id), old_vec.end());
        // 加入到新格子
        cells[new_idx].push_back(site_id);
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

        // --- 修改点：检查是否有足够的列（X, Y, Z, Trace, GradX, GradY, GradZ） ---
        // 现在至少需要7列
        if (values.size() >= 7) { 
            // 提取前四列作为基础数据
            double x = values[0];
            double y = values[1];
            double z = values[2];
            double trace = values[3];

            // --- 修改点：提取后三列作为梯度向量 ---
            // 将后三列数据放入一个 std::vector<double>
            std::vector<double> gradient_vector;
            gradient_vector.push_back(values[4]); // GradX
            gradient_vector.push_back(values[5]); // GradY
            gradient_vector.push_back(values[6]); // GradZ

            // --- 修改点：调用 StressPoint 的新构造函数 ---
            // 传入 gradient_vector
            stress_field_data.emplace_back(current_id++, x, y, z, trace, gradient_vector);
        } else {
            std::cerr << "警告：应力场文件行格式不正确（预期至少7列），跳过此行： " << line << std::endl;
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

int KMC_Simulator::find_closest_Oxy_id(double px, double py, double pz) const {
    // 确定查询点所在的格子坐标
    int cx = static_cast<int>((px) / cell_size_x);
    int cy = static_cast<int>((py) / cell_size_y);
    int cz = static_cast<int>((pz) / cell_size_z);

    double min_dist_sq = std::numeric_limits<double>::max();
    int closest_id = -1;

    // 搜索周围 3x3x3 共 27 个格子
    for (int i = cx - 1; i <= cx + 1; ++i) {
        if (i < 0 || i >= nx) continue;
        for (int j = cy - 1; j <= cy + 1; ++j) {
            if (j < 0 || j >= ny) continue;
            for (int k = cz - 1; k <= cz + 1; ++k) {
                if (k < 0 || k >= nz) continue;

                int cell_idx = i + j * nx + k * nx * ny;
                // 只看这一个格子里的一小撮原子
                for (int sid : cells[cell_idx]) {
                    if (sites[sid].type == OXYG) { // 找到氧原子
                        double dx = px - sites[sid].x;
                        double dy = py - sites[sid].y;
                        double dz = pz - sites[sid].z;
                        double d2 = dx*dx + dy*dy + dz*dz;
                        if (d2 < min_dist_sq) {
                            min_dist_sq = d2;
                            closest_id = sid;
                        }
                    }
                }
            }
        }
    }
    return closest_id;
}

int KMC_Simulator::find_closest_Non_Oxy_id(double px, double py, double pz) const {
    int cx = static_cast<int>((px) / cell_size_x);
    int cy = static_cast<int>((py) / cell_size_y);
    int cz = static_cast<int>((pz) / cell_size_z);

    double min_dist_sq = std::numeric_limits<double>::max();
    int closest_id = -1;

    // 搜索周围 3x3x3 共 27 个格子
    for (int i = cx - 1; i <= cx + 1; ++i) {
        if (i < 0 || i >= nx) continue;
        for (int j = cy - 1; j <= cy + 1; ++j) {
            if (j < 0 || j >= ny) continue;
            for (int k = cz - 1; k <= cz + 1; ++k) {
                if (k < 0 || k >= nz) continue;

                int cell_idx = i + j * nx + k * nx * ny;
                // 只看这一个格子里的一小撮原子
                for (int sid : cells[cell_idx]) {
                    if (sites[sid].type != OXYG) { // 找到氧原子
                        double dx = px - sites[sid].x;
                        double dy = py - sites[sid].y;
                        double dz = pz - sites[sid].z;
                        double d2 = dx*dx + dy*dy + dz*dz;
                        if (d2 < min_dist_sq) {
                            min_dist_sq = d2;
                            closest_id = sid;
                        }
                    }
                }
            }
        }
    }
    return closest_id;
}

double KMC_Simulator::calculate_adsorption_propensity() const {
    const double ADSORPTION_RATE_PER_SECOND = 2.0e3; 

    if (current_time < 1.0) { 
        return 1e-12; // 初始阶段几乎没有氧气吸附
    }
    int oxygen_count = 0;
    for (const auto& s : sites) {
        if (s.type == OXYG) {
            oxygen_count++;
        }
    }
    if (oxygen_count >= 6000) { 
        return 1.1e-12;
    }
    double characteristic_oxithickness = 0.2;
    double total_oxide_thickness = num_ni_oxide * NI_OXIDE_THICKNESS + num_si_oxide * SI_OXIDE_THICKNESS + num_cr_oxide * CR_OXIDE_THICKNESS;
    double oxide_factor =4.0 - 3.0*std::exp(-total_oxide_thickness / characteristic_oxithickness);
    double characteristic_stress = 1000.0;
    double stress_factor = 4.0 - 3.0*std::exp(-GB_normal_stress / characteristic_stress);
    return ADSORPTION_RATE_PER_SECOND*oxide_factor*stress_factor; 
}


void KMC_Simulator::calculate_all_propensities_and_events() {
    total_rate=0.0;
    all_events.clear();                 
    event_propensities_for_selection.clear(); 
    site_to_event_indices.clear();      

    int current_event_index = 0; // 用于跟踪当前事件在 all_events 向量中的索引
    double adsorption_propensity = calculate_adsorption_propensity();
    if (adsorption_propensity > 0.0) {
        // 使用一个哑元 (Dummy) Site 对象
        // static Site dummy_adsorption_site(std::numeric_limits<int>::max(), 0.0, 0.0, 0.0, -1); 

        Event adsorption_event(OXYGEN_ADSORPTION_EVENT, adsorption_propensity, -1);
        all_events.push_back(adsorption_event);
        event_propensities_for_selection.push_back(adsorption_propensity);
        current_event_index++; 
    }
    for (int i = 0; i < num_sites; ++i) {
        // **计算事件倾向：**
        // 对于随机游走事件，我们假设其倾向是一个固定值，例如 1.0。
        // 如果未来引入依赖于 Site 状态（如是否有空位）的倾向，这部分逻辑会更复杂。
        if (sites[i].type==SILIC){
        double calculated_propensity = calculate_Si_site_random_walk_propensity(sites[i]); 
        Event rw_event(SILI_DIFF, calculated_propensity, sites[i].id); 
        // 将创建的 Event 对象添加到总事件列表 `all_events`
        all_events.push_back(rw_event); 
        event_propensities_for_selection.push_back(calculated_propensity); 
        site_to_event_indices[sites[i].id].push_back(current_event_index);
        current_event_index++;}
        else if (sites[i].type==CROM)
        {
        double calculated_propensity = calculate_Cr_site_random_walk_propensity(sites[i]); 
        Event rw_event(CROM_DIFF, calculated_propensity, sites[i].id); 
        // 将创建的 Event 对象添加到总事件列表 `all_events`
        all_events.push_back(rw_event); 
        event_propensities_for_selection.push_back(calculated_propensity); 
        site_to_event_indices[sites[i].id].push_back(current_event_index);
        current_event_index++;
        }
        else if (sites[i].type==NIK)
        {
        double calculated_propensity = calculate_Ni_site_random_walk_propensity(sites[i]); 
        Event rw_event(NIK_DIFF, calculated_propensity, sites[i].id); 
        // 将创建的 Event 对象添加到总事件列表 `all_events`
        all_events.push_back(rw_event); 
        event_propensities_for_selection.push_back(calculated_propensity); 
        site_to_event_indices[sites[i].id].push_back(current_event_index);
        current_event_index++;
        }
        total_rate=std::accumulate(event_propensities_for_selection.begin(), event_propensities_for_selection.end(), 0.0); 
    }
    // rate_scale = std::accumulate(event_propensities_for_selection.begin(), event_propensities_for_selection.end(), 0.0);
}
double KMC_Simulator::calculate_oxy_diffusion_propensity(double x,double y,double z){
    double pre_factor=5e2;
    double charateristic_stress=1500.0;
    double stress_factor=3.0-2.0*std::exp(-GB_normal_stress/charateristic_stress);
    double oxy_diff_propensity=pre_factor*stress_factor*calculate_prop_scalar(x,y,z);
    return oxy_diff_propensity;
}

double KMC_Simulator::calculate_Si_site_random_walk_propensity(const Site& s) const {
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
    // double exponent_term = effective_barrier / kb_T;
    double dis_to_GB=grain_boundary.get_GB_distance(s.x, s.y, s.z);
    double effective_chemical_potential=kb_T*DOSE*0.15*std::log(1/10.0)*std::exp(-dis_to_GB/5e-6);
    double exponent_term = (effective_barrier - effective_chemical_potential)/ kb_T;
    // 最终计算倾向： v0 * exp(-exponent_term)
    double propensity = PRE_FACTOR_V0 * std::exp(-exponent_term);
    propensity *= 1e-3; 
    // ms-1 unit
    propensity *= calculate_prop_scalar(s.x, s.y, s.z);

    // 4. 确保计算出的倾向为正数，防止数学错误或不合理结果
    if (propensity <= 0.0) {
        propensity = 1e-30; // 设置一个非常小的正值，避免零或负倾向
        std::cerr << "警告：为 Site " << s.id << " 计算出的倾向为非正值 (" << propensity << ")，已强制设置为 1e-30。最近应力迹: " << site_stress_trace << std::endl;
    }

    return propensity;
}
double KMC_Simulator::calculate_Cr_site_random_walk_propensity(const Site& s) const {
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
   
    double dis_to_GB=grain_boundary.get_GB_distance(s.x, s.y, s.z);
    double effective_chemical_potential=kb_T*DOSE*0.15*std::log(10.0)*std::exp(-dis_to_GB/5e-6);
    double exponent_term = (effective_barrier - effective_chemical_potential)/ kb_T;
    // 最终计算倾向： v0 * exp(-exponent_term)
    double propensity = PRE_FACTOR_V0 * std::exp(-exponent_term);
    propensity *= 1e-3; 
    // ms-1 unit
    propensity *= calculate_prop_scalar(s.x, s.y, s.z);

    // 4. 确保计算出的倾向为正数，防止数学错误或不合理结果
    if (propensity <= 0.0) {
        propensity = 1e-30; // 设置一个非常小的正值，避免零或负倾向
        std::cerr << "警告：为 Site " << s.id << " 计算出的倾向为非正值 (" << propensity << ")，已强制设置为 1e-30。最近应力迹: " << site_stress_trace << std::endl;
    }

    return propensity;
}

double KMC_Simulator::calculate_Ni_site_random_walk_propensity(const Site& s) const {
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
    // double exponent_term = effective_barrier / kb_T;
    double dis_to_GB=grain_boundary.get_GB_distance(s.x, s.y, s.z);
    double effective_chemical_potential=kb_T*DOSE*0.15*std::log(1/10.0)*std::exp(-dis_to_GB/5e-6);
    double exponent_term = (effective_barrier - effective_chemical_potential)/ kb_T;
    // 最终计算倾向： v0 * exp(-exponent_term)
    double propensity = PRE_FACTOR_V0 * std::exp(-exponent_term);
    propensity *= 1e-3; // ms-1 unit
    propensity *= calculate_prop_scalar(s.x, s.y, s.z);
   
    // 4. 确保计算出的倾向为正数，防止数学错误或不合理结果
    if (propensity <= 0.0) {
        propensity = 1e-30; // 设置一个非常小的正值，避免零或负倾向
        std::cerr << "警告：为 Site " << s.id << " 计算出的倾向为非正值 (" << propensity << ")，已强制设置为 1e-30。最近应力迹: " << site_stress_trace << std::endl;
    }

    return propensity;
}

double KMC_Simulator::calculate_oxi_prob(int type){
    if (current_time<=1.0){
        return 1e-15;
    }
    else {
        double prob;
        double preprob;
        if (type==SILIC){
            preprob=std::exp(-(OXIDATION_BARRIER+ALPHA*SI_OXIDATION_GIBBS)/(KB*TEMPERATURE));
        }
        else if (type==CROM)
        {
            preprob=std::exp(-(OXIDATION_BARRIER+ALPHA*CR_OXIDATION_GIBBS)/(KB*TEMPERATURE));
        }
        else if (type==NIK)
        {
            preprob=std::exp(-(OXIDATION_BARRIER+ALPHA*NI_OXIDATION_GIBBS)/(KB*TEMPERATURE));
        }
        else {
            preprob=1e-15;
        }
        prob=preprob/(1.0+std::exp((num_cr_oxide-100.0)/3.0));
        double p = 1.0-std::pow(1-prob,rate_scale/total_rate);
        return p;
    }
}

double KMC_Simulator::calculate_jump_scalar(double x,double y,double z) const{
    double dis_to_GB=grain_boundary.get_GB_distance(x,y,z);
    double scalar;
    if (dis_to_GB<3e-5){
        scalar=1.0*std::pow(10.0,(1.0-2.0)*(1.0/(1.0+std::exp((dis_to_GB-3e-5)/5e-6))));
    }
    else {
        scalar=1.0*std::pow(10.0,(1.0-2.0)*(1.0/(1.0+std::exp((dis_to_GB-3e-5)/1e-4))));
    }
    return scalar;
}

double KMC_Simulator::calculate_prop_scalar(double x,double y,double z) const{
    double scalar=calculate_jump_scalar(x,y,z);
    double p=1.0/(scalar*scalar);
    return p;
}

double KMC_Simulator::normal_sample(double mean, double stddev)
{
    static thread_local std::mt19937_64 rng{std::random_device{}()};
    std::normal_distribution<double> dist(mean, stddev);
    return dist(rng);
    // return mean;
}

int KMC_Simulator::select_event_index() {
    // 计算所有事件的总倾向 H (吉莱斯皮算法中的 A_0)
    // double total_propensity = std::accumulate(event_propensities_for_selection.begin(), event_propensities_for_selection.end(), 0.0);
    double total_propensity = total_rate;
    // 如果总倾向为零或负，说明没有事件可以发生
    if (total_propensity <= 0.0) {
        std::cout<<"total p < 0"<<std::endl;
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
    

    // Site& target_site = const_cast<Site&>(sites[chosen_event.site_id]); 

    // 根据 Event 对象的类型 (etype) 来执行不同的操作
    // KMC_Simulator 现在承担了所有事件类型的执行逻辑。
    switch (chosen_event.etype) {
        case SILI_DIFF: {
            Site& target_site = const_cast<Site&>(sites[chosen_event.site_id]); 
            int closest_sp_id = find_closest_stress_point_id(target_site.x, target_site.y, target_site.z);
            double ox=target_site.x, oy=target_site.y, oz=target_site.z;
            std::vector<double> direction(3);
            double displacement;
            if (Diff_mod==TRACKDIFF){
                direction=stress_field_data[closest_sp_id].tr_grad;
                double displacement = jump_distance;
            }
            else if (Diff_mod==RANDDIFF)
            {
                double adjusted_jump_distance = jump_distance*calculate_jump_scalar(target_site.x,target_site.y,target_site.z);
                double rx = normal_sample(0.0,adjusted_jump_distance);
                double ry = normal_sample(0.0,adjusted_jump_distance);
                double rz = normal_sample(0.0,adjusted_jump_distance);
                double leng = std::sqrt(rx*rx+ry*ry+rz*rz);
                direction[0]=rx/leng;
                direction[1]=ry/leng;
                direction[2]=rz/leng;
                displacement = leng;
            }
            

            // 调用 Site 对象的 `move` 方法来更新其位置
            target_site.move(direction, displacement);
            apply_pbc(target_site); 
            update_site_cell(target_site.id, ox, oy, oz);
            total_rate-=all_events[event_index].propensity;
            double updated_prop=calculate_Si_site_random_walk_propensity(target_site);
            all_events[event_index].propensity=updated_prop;
            event_propensities_for_selection[event_index]=updated_prop;
            total_rate+=updated_prop;
            int Oxyid = find_closest_Oxy_id(target_site.x, target_site.y, target_site.z);
            if (Oxyid!=-1){
                double Oxy_dis=std::sqrt((target_site.x-sites[Oxyid].x)*(target_site.x-sites[Oxyid].x)+(target_site.y-sites[Oxyid].y)*(target_site.y-sites[Oxyid].y)+(target_site.z-sites[Oxyid].z)*(target_site.z-sites[Oxyid].z));
                if (Oxy_dis<OXIDATION_THRESH && sites[Oxyid].status==0 && target_site.status==0 && get_uniform_random()<calculate_oxi_prob(target_site.type)){
                    updated_prop=1e-15;
                    total_rate-=all_events[event_index].propensity;
                    total_rate+=updated_prop;
                    all_events[event_index].propensity=updated_prop;
                    event_propensities_for_selection[event_index]=updated_prop;
                    
                    int Oxy_diff_index=-1;
                    for (int i=0; i<site_to_event_indices[Oxyid].size(); ++i){
                        if (all_events[site_to_event_indices[Oxyid][i]].etype==OXYGEN_DIFF){
                            Oxy_diff_index=site_to_event_indices[Oxyid][i];
                            break;
                        }

                    }
                    total_rate-=all_events[Oxy_diff_index].propensity;
                    total_rate+=updated_prop;
                    all_events[Oxy_diff_index].propensity=updated_prop;
                    event_propensities_for_selection[Oxy_diff_index]=updated_prop;
                    std::cout<<"Si oxide formed"<<std::endl;
                    num_si_oxide += 1;
                    num_total_oxide += 1;
                    target_site.status=1;
                    sites[Oxyid].status=1;
                    write_oxide_csv();

                }
            }
            break;
        }
        case CROM_DIFF: {
            Site& target_site = const_cast<Site&>(sites[chosen_event.site_id]); 
            int closest_sp_id = find_closest_stress_point_id(target_site.x, target_site.y, target_site.z);
            double ox=target_site.x, oy=target_site.y, oz=target_site.z;
            std::vector<double> direction(3);
            double displacement;
            if (Diff_mod==TRACKDIFF){
                direction=stress_field_data[closest_sp_id].tr_grad;
                double displacement = jump_distance;
            }
            else if (Diff_mod==RANDDIFF)
            {
                double adjusted_jump_distance = jump_distance*calculate_jump_scalar(target_site.x,target_site.y,target_site.z);
                double rx = normal_sample(0.0,adjusted_jump_distance);
                double ry = normal_sample(0.0,adjusted_jump_distance);
                double rz = normal_sample(0.0,adjusted_jump_distance);
                double leng = std::sqrt(rx*rx+ry*ry+rz*rz);
                direction[0]=rx/leng;
                direction[1]=ry/leng;
                direction[2]=rz/leng;
                displacement = leng;
            }
            // 调用 Site 对象的 `move` 方法来更新其位置
            target_site.move(direction, displacement);
            apply_pbc(target_site); 
            update_site_cell(target_site.id, ox, oy, oz);
            double updated_prop=calculate_Cr_site_random_walk_propensity(target_site);
            total_rate-=all_events[event_index].propensity;
                    total_rate+=updated_prop;
            all_events[event_index].propensity=updated_prop;
            event_propensities_for_selection[event_index]=updated_prop;
            int Oxyid = find_closest_Oxy_id(target_site.x, target_site.y, target_site.z);
            if (Oxyid!=-1){
                double Oxy_dis=std::sqrt((target_site.x-sites[Oxyid].x)*(target_site.x-sites[Oxyid].x)+(target_site.y-sites[Oxyid].y)*(target_site.y-sites[Oxyid].y)+(target_site.z-sites[Oxyid].z)*(target_site.z-sites[Oxyid].z));
                if (Oxy_dis<OXIDATION_THRESH && sites[Oxyid].status==0 && target_site.status==0 && get_uniform_random()<calculate_oxi_prob(target_site.type)){
                    updated_prop=1e-15;
                    total_rate-=all_events[event_index].propensity;
                    total_rate+=updated_prop;
                    all_events[event_index].propensity=updated_prop;
                    event_propensities_for_selection[event_index]=updated_prop;
                    int Oxy_diff_index=-1;
                    for (int i=0; i<site_to_event_indices[Oxyid].size(); ++i){
                        if (all_events[site_to_event_indices[Oxyid][i]].etype==OXYGEN_DIFF){
                            Oxy_diff_index=site_to_event_indices[Oxyid][i];
                            break;
                        }

                    }
                    total_rate-=all_events[Oxy_diff_index].propensity;
                    total_rate+=updated_prop;
                    all_events[Oxy_diff_index].propensity=updated_prop;
                    event_propensities_for_selection[Oxy_diff_index]=updated_prop;
                    std::cout<<"Cr oxide formed"<<std::endl;
                    num_cr_oxide += 1;
                    num_total_oxide += 1;
                    target_site.status=1;
                    sites[Oxyid].status=1;
                    write_oxide_csv();
                    

                }
            }
            break;
        }
        
        case NIK_DIFF: {
            Site& target_site = const_cast<Site&>(sites[chosen_event.site_id]); 
            int closest_sp_id = find_closest_stress_point_id(target_site.x, target_site.y, target_site.z);
            double ox=target_site.x, oy=target_site.y, oz=target_site.z;
            std::vector<double> direction(3);
            double displacement;
            if (Diff_mod==TRACKDIFF){
                direction=stress_field_data[closest_sp_id].tr_grad;
                double displacement = jump_distance;
            }
            else if (Diff_mod==RANDDIFF)
            {
                double adjusted_jump_distance = jump_distance*calculate_jump_scalar(target_site.x,target_site.y,target_site.z);
                double rx = normal_sample(0.0,adjusted_jump_distance);
                double ry = normal_sample(0.0,adjusted_jump_distance);
                double rz = normal_sample(0.0,adjusted_jump_distance);
                double leng = std::sqrt(rx*rx+ry*ry+rz*rz);
                direction[0]=rx/leng;
                direction[1]=ry/leng;
                direction[2]=rz/leng;
                displacement = leng;
            }
            // 调用 Site 对象的 `move` 方法来更新其位置
            target_site.move(direction, displacement);
            apply_pbc(target_site); 
            update_site_cell(target_site.id, ox, oy, oz);
            
            double updated_prop=calculate_Ni_site_random_walk_propensity(target_site);
            total_rate-=all_events[event_index].propensity;
            total_rate+=updated_prop;
            all_events[event_index].propensity=updated_prop;
            event_propensities_for_selection[event_index]=updated_prop;
            int Oxyid = find_closest_Oxy_id(target_site.x, target_site.y, target_site.z);
            if (Oxyid!=-1){
                double Oxy_dis=std::sqrt((target_site.x-sites[Oxyid].x)*(target_site.x-sites[Oxyid].x)+(target_site.y-sites[Oxyid].y)*(target_site.y-sites[Oxyid].y)+(target_site.z-sites[Oxyid].z)*(target_site.z-sites[Oxyid].z));
                if (Oxy_dis<OXIDATION_THRESH && sites[Oxyid].status==0 && target_site.status==0 && get_uniform_random()<calculate_oxi_prob(target_site.type)){
                    updated_prop=1e-15;
                    total_rate-=all_events[event_index].propensity;
                    total_rate+=updated_prop;
                    all_events[event_index].propensity=updated_prop;
                    event_propensities_for_selection[event_index]=updated_prop;
                    int Oxy_diff_index=-1;
                    for (int i=0; i<site_to_event_indices[Oxyid].size(); ++i){
                        if (all_events[site_to_event_indices[Oxyid][i]].etype==OXYGEN_DIFF){
                            Oxy_diff_index=site_to_event_indices[Oxyid][i];
                            break;
                        }

                    }
                    total_rate-=all_events[Oxy_diff_index].propensity;
                    total_rate+=updated_prop;
                    all_events[Oxy_diff_index].propensity=updated_prop;
                    event_propensities_for_selection[Oxy_diff_index]=updated_prop;
                    std::cout<<"Ni oxide formed"<<std::endl;
                    num_ni_oxide += 1;
                    num_total_oxide += 1;
                    target_site.status=1;
                    sites[Oxyid].status=1;
                    write_oxide_csv();

                }
            }
            break;
        }

        case OXYGEN_ABS: {
            int new_site_id = 0;
            if (!sites.empty()) {
                // 查找当前 sites 列表中最大的 ID + 1
                int max_id = sites[0].id;
                for(const auto& s : sites) {
                    if (s.id > max_id) max_id = s.id;
                }
                new_site_id = max_id + 1;
            }
            
            // 使用 GB 类在交线上随机生成 OXYGEN Site
            Site new_oxygen = grain_boundary.generate_random_oxygen_site(new_site_id);
            
            // 将新 Site 加入到 sites 列表中
            sites.push_back(new_oxygen);
            Event oxy_diff_event(OXYGEN_DIFF, calculate_oxy_diffusion_propensity(new_oxygen.x,new_oxygen.y,new_oxygen.z),new_oxygen.id);
            all_events.push_back(oxy_diff_event);
            event_propensities_for_selection.push_back(calculate_oxy_diffusion_propensity(new_oxygen.x,new_oxygen.y,new_oxygen.z));
            total_rate+=calculate_oxy_diffusion_propensity(new_oxygen.x,new_oxygen.y,new_oxygen.z);
            // total_rate-=all_events[event_index].propensity;
            // total_rate+=calculate_adsorption_propensity();
            // all_events[event_index].propensity=calculate_adsorption_propensity();
            // event_propensities_for_selection[event_index]=all_events[event_index].propensity;
            int current_event_index=all_events.size()-1;
            site_to_event_indices[new_oxygen.id].push_back(current_event_index);
            std::cout << "Adsorption: New Oxygen Site " << new_site_id 
                      << " created at (" << new_oxygen.x << ", " << new_oxygen.y 
                      << ", " << new_oxygen.z << ") on GB." << std::endl;
            double absoption_propensity=calculate_adsorption_propensity();
            total_rate-=all_events[event_index].propensity;
            total_rate+=absoption_propensity;
            all_events[event_index].propensity=absoption_propensity;
            event_propensities_for_selection[event_index]=all_events[event_index].propensity; 
            break;

        }
        case OXYGEN_DIFF: {
            Site& target_site = const_cast<Site&>(sites[chosen_event.site_id]);
            double ox=target_site.x, oy=target_site.y, oz=target_site.z;
            double xm = 2.0*get_uniform_random()-1;
            double zm = 2.0*get_uniform_random()-1;
            target_site.x += 10.0*xm*jump_distance*calculate_jump_scalar(target_site.x,target_site.y,target_site.z);
            target_site.z += 10.0*zm*jump_distance*calculate_jump_scalar(target_site.x,target_site.y,target_site.z);
            target_site.y = -(GB_A*target_site.x+GB_C*target_site.z+GB_D)/GB_B;
            grain_boundary.apply_pbc_and_constrain(target_site);
            update_site_cell(target_site.id, ox, oy, oz);


            // std::cout<<"oxygen move, site id: "<<target_site.id<<" "<<jump_distance<<std::endl;
            int NonOxyid = find_closest_Non_Oxy_id(target_site.x, target_site.y, target_site.z);
            if (NonOxyid!=-1){
                double Oxy_dis=std::sqrt((target_site.x-sites[NonOxyid].x)*(target_site.x-sites[NonOxyid].x)+(target_site.y-sites[NonOxyid].y)*(target_site.y-sites[NonOxyid].y)+(target_site.z-sites[NonOxyid].z)*(target_site.z-sites[NonOxyid].z));
                if (Oxy_dis<OXIDATION_THRESH && sites[NonOxyid].status==0 && target_site.status==0 && get_uniform_random()<calculate_oxi_prob(sites[NonOxyid].type)){
                    double updated_prop=1e-15;
                    
                    total_rate-=all_events[event_index].propensity;
                    total_rate+=updated_prop;
                    all_events[event_index].propensity=updated_prop;
                    event_propensities_for_selection[event_index]=updated_prop;
                    int NonOxy_diff_index=-1;
                    for (int i=0; i<site_to_event_indices[NonOxyid].size(); ++i){
                        if (all_events[site_to_event_indices[NonOxyid][i]].etype==SILI_DIFF || all_events[site_to_event_indices[NonOxyid][i]].etype==CROM_DIFF || all_events[site_to_event_indices[NonOxyid][i]].etype==NIK_DIFF){
                            NonOxy_diff_index=site_to_event_indices[NonOxyid][i];
                            break;
                        }

                    }
                    total_rate-=all_events[NonOxy_diff_index].propensity;
                    total_rate+=updated_prop;
                    all_events[NonOxy_diff_index].propensity=updated_prop;
                    event_propensities_for_selection[NonOxy_diff_index]=updated_prop;
                    if (sites[NonOxyid].type==SILIC){
                        std::cout<<"Si oxide formed"<<std::endl;
                        num_si_oxide += 1;
                        num_total_oxide += 1;
                    }
                    else if (sites[NonOxyid].type==CROM){
                        std::cout<<"Cr oxide formed"<<std::endl;
                        num_cr_oxide += 1;
                        num_total_oxide += 1;
                        // oxi_prob=calculate_oxi_prob();
                    }
                    else if (sites[NonOxyid].type==NIK){
                        std::cout<<"Ni oxide formed"<<std::endl;
                        num_ni_oxide += 1;
                        num_total_oxide += 1;
                    }
                    target_site.status=1;
                    sites[NonOxyid].status=1;
                    write_oxide_csv();

                }
            }
            
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


void KMC_Simulator::print_status() {
    std::cout << "Time: " << std::fixed << std::setprecision(6) << current_time 
              << " | Steps: " << total_steps <<" total rate: "<<total_rate <<std::endl;
}

// --- 输出 Site 坐标到文件 ---
// 将所有 Site 的当前三维坐标以 .xyz 格式写入文件，方便外部工具（如 VMD）进行可视化。
void KMC_Simulator::dump_sites(long long int step) {
    std::string filename = output_dir+"/"+"sites_" + std::to_string(step) + ".xyz"; // 文件名包含 KMC 步数
    std::ofstream outfile(filename); // 创建并打开文件

    if (!outfile.is_open()) { // 检查文件是否成功打开
        std::cerr << "Error: Could not open file " << filename << " for dumping sites." << std::endl;
        return;
    }

    // .xyz 文件格式要求：
    outfile << sites.size() << std::endl; // 第一行：原子数量
    outfile << "KMC Step: " << step << ", Time: " << std::fixed << std::setprecision(6) << current_time << std::endl; // 第二行：注释

    // 写入每个 Site 的数据
    for (const auto& s : sites) {
        if (s.type==SILIC){
        outfile << "Si " << std::fixed << std::setprecision(6) 
                << s.x << " " << s.y << " " << s.z << std::endl;
        }
        else if (s.type==CROM){
        outfile << "Cr " << std::fixed << std::setprecision(6) 
                << s.x << " " << s.y << " " << s.z << std::endl;
        }
        else if (s.type==OXYG){
        outfile << "O " << std::fixed << std::setprecision(6) 
                << s.x << " " << s.y << " " << s.z << std::endl;}
        else if (s.type==NIK){
        outfile << "Ni " << std::fixed << std::setprecision(6) 
                << s.x << " " << s.y << " " << s.z << std::endl;}
        
        // std::cout<<"writing"<<std::endl;
    }

    outfile.close(); // 关闭文件
}

void KMC_Simulator::write_sites_csv() {
    std::string filename = output_dir+"/"+"sitesresult" + output_dir + ".csv"; // 文件名包含 KMC 步数
    std::ofstream outfile(filename); // 创建并打开文件

    if (!outfile.is_open()) { // 检查文件是否成功打开
        std::cerr << "Error: Could not open file " << filename << " for dumping sites." << std::endl;
        return;
    }

    // // .xyz 文件格式要求：
    // outfile << sites.size() << std::endl; // 第一行：原子数量
    // outfile << "KMC Step: " << step << ", Time: " << std::fixed << std::setprecision(6) << current_time << std::endl; // 第二行：注释

    // 写入每个 Site 的数据
    for (const auto& s : sites) {
        if (s.status==0){
        if (s.type==SILIC){
        outfile << "Si," << std::fixed << std::setprecision(6) 
                << s.x << "," << s.y << "," << s.z << std::endl;
        }
        else if (s.type==CROM){
        outfile << "Cr," << std::fixed << std::setprecision(6) 
                << s.x << "," << s.y << "," << s.z << std::endl;
        }
        else if (s.type==OXYG){
        outfile << "O," << std::fixed << std::setprecision(6) 
                << s.x << "," << s.y << "," << s.z << std::endl;}
        else if (s.type==NIK){
        outfile << "Ni," << std::fixed << std::setprecision(6) 
                << s.x << "," << s.y << "," << s.z << std::endl;}
        
        }// std::cout<<"writing"<<std::endl;
    }

    outfile.close(); // 关闭文件
}


void KMC_Simulator::write_oxide_csv(){
    const std::string filename = output_dir+"/"+"oxidegrow"+output_dir+".csv";

    // 以追加模式打开文件（不存在则创建）
    std::ofstream fout(filename, std::ios::app);
    if (!fout.is_open()) {
        std::cerr << "Error: Cannot open " << filename << " for writing.\n";
        return;
    }

    // 关键：判断当前写位置是否在文件开头
    // 如果文件是刚创建的，tellp() 会返回 0，需要写表头
    bool need_header = (fout.tellp() == 0);

    if (need_header) {
        fout << "time,Si_oxide,Cr_oxide,Ni_oxide,Total_oxide\n";
    }

    // 写入一行数据
    fout << current_time << ","
         << num_si_oxide << ","
         << num_cr_oxide << ","
         << num_ni_oxide << ","
         << num_total_oxide << "\n";



}

void KMC_Simulator::write_propensity_csv(){
    const std::string filename = output_dir+"/"+"propensity"+std::to_string(current_time)+".csv";

    // 以追加模式打开文件（不存在则创建）
    std::ofstream fout(filename, std::ios::app);
    if (!fout.is_open()) {
        std::cerr << "Error: Cannot open " << filename << " for writing.\n";
        return;
    }

    // 关键：判断当前写位置是否在文件开头
    // 如果文件是刚创建的，tellp() 会返回 0，需要写表头
    bool need_header = (fout.tellp() == 0);

    if (need_header) {
        fout << "id,type,propensity,propensity_select\n";
    }
    int num_events=event_propensities_for_selection.size();
    for (int i=0;i<num_events;++i){
        fout << i << ","
        <<all_events[i].etype<<","
        <<all_events[i].propensity<<","
        <<event_propensities_for_selection[i]<<"\n";
    }

}

void KMC_Simulator::run(double max_time, long long int max_steps) {
    std::cout << "\nStarting KMC simulation..." << std::endl;
    std::cout << "Maximum Simulation Time: " << max_time 
              << ", Maximum KMC Steps: " << max_steps << std::endl;

    // **首次构建事件列表和计算所有倾向**
    fs::create_directories(output_dir);  // 目录不存在就创建（已存在也不会报错）
    init_cell_list( nx, ny, nz);
    calculate_all_propensities_and_events();
    rate_scale = 1e9;
    // rate_scale = std::accumulate(event_propensities_for_selection.begin(), event_propensities_for_selection.end(), 0.0);
    dump_sites(0); // 输出初始构型
    write_propensity_csv();
    write_oxide_csv();

    // KMC 模拟主循环
    while (current_time < max_time && total_steps < max_steps) {
        // 获取当前所有事件的总倾向 H
        // double total_propensity = std::accumulate(event_propensities_for_selection.begin(), event_propensities_for_selection.end(), 0.0);
        double total_propensity = total_rate;
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
        jump_distance=unit_jump_distance;
        execute_event(chosen_event_idx); 
        total_steps++; // 增加已完成的 KMC 步数
        if (all_events[0].etype==OXYGEN_ABS){
            double old_absorption_propensity=all_events[0].propensity;
            double new_absorption_propensity=calculate_adsorption_propensity();
            total_rate-=old_absorption_propensity;
            total_rate+=new_absorption_propensity;
            event_propensities_for_selection[0]=new_absorption_propensity;
            all_events[0].propensity=new_absorption_propensity;
        }
        else {
            std::cerr << "Warning: First event in all_events is not OXYGEN_ABS. Check event list construction." << std::endl;
            
        }

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
        if (total_steps % 1000000 == 0) { // 例如，每 1000 步输出一次
            print_status();
            dump_sites(total_steps);
            write_sites_csv();
        }
        if (total_steps % 100000000 == 0) { // 例如，每 1000 步输出一次
           
            write_propensity_csv();
        }
        if (total_steps == 1000000) { 
           
            write_propensity_csv();
        }
    }
    // 模拟结束时，输出最终状态和构型
    std::cout << "\nSimulation finished. Final status:" << std::endl;
    print_status();
    dump_sites(total_steps); 
    write_sites_csv();
    write_oxide_csv();
}