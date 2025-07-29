#ifndef STRESS_POINT_H
#define STRESS_POINT_H

#include <cmath> // 用于可能的距离计算

class StressPoint {
public:
    
    int id;
    double x, y, z;
    double trace; // 应力迹

    // 构造函数
    StressPoint(int id, double sx = 0.0, double sy = 0.0, double sz = 0.0, double s_trace = 0.0);

    // 计算到另一个点的欧几里得距离
    double get_distance(double px, double py, double pz) const;
};

#endif