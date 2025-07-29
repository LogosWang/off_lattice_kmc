#include "StressPoint.h"

// 构造函数实现，更新以包含 id
StressPoint::StressPoint(int sid, double sx, double sy, double sz, double s_trace)
    : id(sid), x(sx), y(sy), z(sz), trace(s_trace) {
    // 构造函数体为空，因为成员变量已经在初始化列表中初始化
}

// 计算到另一个点的欧几里得距离的实现保持不变
double StressPoint::get_distance(double px, double py, double pz) const {
    double dx = x - px;
    double dy = y - py;
    double dz = z - pz;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}