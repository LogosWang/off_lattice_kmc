#ifndef GB_H
#define GB_H

#include "Site.h"     // 需要 Site 类的定义
#include <vector>     
#include <cmath>      
#include <iostream>   

// 定义与晶界相关的 Site 类型和 Event 类型，以便在 KMC_Simulator 中使用
const int OXYGEN_SITE_TYPE = 2;      // 假设 type = 1 代表 OXYGEN_SITE
const int OXYGEN_ADSORPTION_EVENT = 2; // 假设 etype = 2 代表 OXYGEN_ADSORPTION
const int OXYGEN_DIFFUSION_EVENT = 3;  // 假设 etype = 3 代表 OXYGEN_DIFFUSION

class GB {
public:
    // 晶界平面方程: A*x + B*y + C*z + D = 0
    double A, B, C, D;
    double box_size;            // 模拟盒子的边长，用于 PBC
    double normal_mag_sq;       // 晶界平面法向量的模的平方
    
    // 交线参数：与 Y 顶部 (y = box_size) 的交线段的两个端点
    std::vector<double> intersection_point1; // (x1, box_size, z1)
    std::vector<double> intersection_point2; // (x2, box_size, z2)


    // --- 构造函数 ---
    GB(double a, double b, double c, double d, double b_size);

    // --- 初始化/设置函数 ---

    // 找出 GB 平面与 Y 顶部 (y = box_size) 的交线段端点
    void calculate_z_top_intersection();

    // --- 事件相关函数 ---

    // 在晶界与 Y 顶部交线上的随机位置，生成一个 OXYGEN Site 的初始位置
    Site generate_random_oxygen_site(int new_site_id) const;

    // --- 运动和 PBC 函数 ---

    // 针对 OXYGEN_SITE 的特殊 PBC (迭代求解直到在 Box 内)
    void apply_pbc_and_constrain(Site& s);

    // 将一个位移向量投影到 GB 平面上 (用于 OXYGEN 扩散的位移向量)
    std::vector<double> project_displacement(const std::vector<double>& displacement) const;

private:
    // --- 辅助函数：根据平面方程计算第三个坐标 (用于迭代 PBC) ---

    // 根据 y, z 计算 x，要求 |A| != 0
    double calculate_x(double y, double z) const;
    
    // 根据 x, z 计算 y，要求 |B| != 0
    double calculate_y(double x, double z) const;
    
    // 根据 x, y 计算 z，要求 |C| != 0
    double calculate_z(double x, double y) const;
};

#endif // GB_H