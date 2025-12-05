#include "GB.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

// --- 构造函数 ---
GB::GB(double a, double b, double c, double d, double b_size)
    : A(a), B(b), C(c), D(d), box_size(b_size), intersection_point1(3), intersection_point2(3) {
    
    // 计算法向量模的平方，用于投影（如果需要）
    normal_mag_sq = A*A + B*B + C*C;
    
    if (normal_mag_sq < 1e-12) {
        std::cerr << "Error: GB normal vector magnitude is zero. This is not a valid plane." << std::endl;
    }

    // 初始化交线端点
    calculate_z_top_intersection(); 
}

// --- 私有计算函数实现 ---

double GB::calculate_x(double y, double z) const {
    if (std::abs(A) < 1e-12) {
        // 晶界平行于 X 轴，不能通过 y, z 唯一确定 x
        std::cerr << "Error: GB::calculate_x called with A near zero. Plane is parallel to X-axis." << std::endl;
        return 0.0; 
    }
    // x = -(B*y + C*z + D) / A
    return -(B * y + C * z + D) / A;
}

double GB::calculate_y(double x, double z) const {
    if (std::abs(B) < 1e-12) {
        // 晶界平行于 Y 轴，不能通过 x, z 唯一确定 y
        std::cerr << "Error: GB::calculate_y called with B near zero. Plane is parallel to Y-axis." << std::endl;
        return 0.0; 
    }
    // y = -(A*x + C*z + D) / B
    return -(A * x + C * z + D) / B;
}

double GB::calculate_z(double x, double y) const {
    if (std::abs(C) < 1e-12) {
        // 晶界平行于 Z 轴，不能通过 x, y 唯一确定 z
        std::cerr << "Error: GB::calculate_z called with C near zero. Plane is parallel to Z-axis." << std::endl;
        return 0.0; 
    }
    // z = -(A*x + B*y + D) / C
    return -(A * x + B * y + D) / C;
}


// --- 找出 GB 平面与 Y 顶部 (y = box_size) 的交线段端点 ---
void GB::calculate_z_top_intersection() {
    // 简化处理: 找到与 x=0 和 x=box_size 边界的交点，这两个点定义了交线段
    
    // 1. 尝试找到 x=0 时的交点
    // double z_at_x0 = calculate_z(0.0, box_size);
    // if (z_at_x0 >= 0 && z_at_x0 <= box_size) {
    //     intersection_point1 = {0.0, box_size, z_at_x0};
    // } else {
    //     // 如果 x=0 时的 z 超出边界，使用 z=0 时的交点
    //     double x_at_z0 = calculate_x(box_size, 0.0);
    //     intersection_point1 = {x_at_z0, box_size, 0.0};
    // }
    
    // // 2. 尝试找到 x=box_size 时的交点
    // double z_at_x_max = calculate_z(box_size, box_size);
    // if (z_at_x_max >= 0 && z_at_x_max <= box_size) {
    //     intersection_point2 = {box_size, box_size, z_at_x_max};
    // } else {
    //     // 如果 x=box_size 时的 z 超出边界，使用 z=box_size 时的交点
    //     double x_at_z_max = calculate_x(box_size, box_size);
    //     intersection_point2 = {x_at_z_max, box_size, box_size};
    // }
    
    // 注意：这里的逻辑假设平面与 Y 顶部有一个有效的、跨越 box 边界的交线段。
    // 如果平面完全不接触 y=box_size 表面（例如 y < box_size），这两个点将是无效的，但
    // 对于 KMC 模拟器，我们假定晶界是贯穿的。
    double lineA=A;
    double lineB=B;
    double lineC=C*box_size+D;
    // 交线方程lineAX+lineBY+lineC=0
    double x1 = 0.0;
    double x2 = calculate_x(0.0,box_size);
    double x3 = calculate_x(box_size,box_size);
    double x4 = box_size;
    std::cout<<x1<<" "<<x2<<" "<<x3<<" "<<x4<<std::endl;
    std::cout<<lineA*lineB<<std::endl;
    if (lineA*lineB<0){
        intersection_point1[0]=std::max(x1,x2);
        intersection_point1[2]=box_size;
        intersection_point1[1]=calculate_y(intersection_point1[0],intersection_point1[2]);

        intersection_point2[0]=std::min(x3,x4);
        intersection_point2[2]=box_size;
        intersection_point2[1]=calculate_y(intersection_point2[0],intersection_point2[2]);
        
    }
    else{
        intersection_point1[0]=std::max(x1,x3);
        intersection_point1[2]=box_size;
        intersection_point1[1]=calculate_y(intersection_point1[0],intersection_point1[2]);

        intersection_point2[0]=std::min(x2,x4);
        intersection_point2[2]=box_size;
        intersection_point2[1]=calculate_y(intersection_point2[0],intersection_point2[2]);
        
    }
    std::cout<<"absorption range"<<intersection_point1[0]
    <<" "<<intersection_point1[1]<<" "<<intersection_point1[2]
    <<" to "<<intersection_point2[0]
    <<" "<<intersection_point2[1]<<" "<<intersection_point2[2]
    <<std::endl;
}


// --- 在交线上随机生成一个 OXYGEN Site 的初始位置 ---
Site GB::generate_random_oxygen_site(int new_site_id) const {
    // 1. 计算交线段的向量 V = P2 - P1
    double dx = intersection_point2[0] - intersection_point1[0];
    double dy = intersection_point2[1] - intersection_point1[1];
    double dz = intersection_point2[2] - intersection_point1[2];

    // 2. 生成一个随机参数 t in [0, 1]
    double t = static_cast<double>(rand()) / RAND_MAX;

    // 3. 计算随机点 P_rand = P1 + t * V
    double x = intersection_point1[0] + t * dx;
    double y = intersection_point1[1] + t * dy;
    double z = intersection_point1[2] + t * dz;
    
    // // 确保 y = box_size (在 calculate_y_top_intersection 中已设置)
    // y = box_size; 

    return Site(new_site_id, x, y, z, 2);
}

// --- 将一个位移向量投影到 GB 平面上 ---
std::vector<double> GB::project_displacement(const std::vector<double>& displacement) const {
    // V_parallel = V - (V . N) * N / |N|^2
    
    if (normal_mag_sq < 1e-12) {
        return displacement; // 异常情况，返回原始位移
    }

    // 1. 计算 V . N (点积)
    double dot_product = displacement[0] * A + displacement[1] * B + displacement[2] * C;

    // 2. 计算 V_normal (法向分量)
    double scale = dot_product / normal_mag_sq;

    std::vector<double> v_parallel(3);
    // 3. 计算 V_parallel = V - V_normal
    v_parallel[0] = displacement[0] - scale * A;
    v_parallel[1] = displacement[1] - scale * B;
    v_parallel[2] = displacement[2] - scale * C;

    return v_parallel;
}

double GB::get_GB_distance(double x, double y, double z) const {
    double d=std::abs(A*x+B*y+C*z+D)/std::sqrt(A*A + B*B + C*C);
    return d;
}

// --- 针对 OXYGEN_SITE 的特殊 PBC 和平面约束 (迭代直到收敛) ---
void GB::apply_pbc_and_constrain(Site& s) {
    if (s.type != OXYGEN_SITE_TYPE) {
        return; 
    }

    // 您的要求：无限循环，直到所有坐标都在边界内
    int iter=0;
    while (true) {
        bool boundary_crossed = false;
        
        // --- 1. 检查并调整 X 坐标 ---
        if (s.x < 0 || s.x > box_size) {
            // PBC 调整 X
            if (s.x < 0) s.x += box_size;
            s.x = fmod(s.x, box_size); 
            // if (s.x < 0) s.x += box_size;
            
            // 保持 Y 不变，根据平面方程重新计算 Z
            s.z = calculate_z(s.x, s.y);
            boundary_crossed = true;
        }

        // --- 2. 检查并调整 Y 坐标 ---
        if (s.y < 0 || s.y > box_size) {
            // PBC 调整 Y
            if (s.y < 0) s.y += box_size;
            s.y = fmod(s.y, box_size); 
            // if (s.y < 0) s.y += box_size;
            
            // 保持 Z 不变，根据平面方程重新计算 X
            s.x = calculate_x(s.y, s.z);
            boundary_crossed = true;
        }

        // --- 3. 检查并调整 Z 坐标 ---
        if (s.z < 0 || s.z > box_size) {
            // PBC 调整 Z
            if (s.z < 0) s.z += box_size;
            s.z = fmod(s.z, box_size); 
            // if (s.z < 0) s.z += box_size;
            
            // 保持 X 不变，根据平面方程重新计算 Y
            s.y = calculate_y(s.x, s.z);
            boundary_crossed = true;
        }

        // 如果本轮循环中没有发生任何边界调整，说明点已经同时在盒子内和平面上，退出循环
        iter ++;
        if (iter>50){
            std::cout<<"2D PBC ecceed maximum iteration number."<<std::endl;
            break;
        }
        if (!boundary_crossed) {
            // std::cout<<"iteration number:"<<iter<<std::endl;
            break;
        }
    }
}