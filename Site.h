#ifndef SITE_H
#define SITE_H

// 将struct改为class
class Site {
public: // 通常情况下，站点的属性（位置、ID、类型）会设置为public以便外部直接访问，特别是对于简单的数据结构
    int id;       // 站点ID
    double x, y, z; // 坐标
    int type;     // 站点类型 (如果需要)

    // 构造函数
    Site(int site_id, double initial_x, double initial_y, double initial_z, int site_type = 0);
    void move(int,double);
    // 示例：可以添加一些成员函数
    // double get_distance(const Site& other) const; // 计算与另一个站点的距离
    // void move_by(double dx, double dy, double dz); // 根据位移向量移动站点
};

#endif // SITE_H