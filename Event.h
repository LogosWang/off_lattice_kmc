#ifndef EVENT_H
#define EVENT_H

// 包含 Site 类的定义，因为事件的执行和倾向计算可能需要 Site 对象
#include "Site.h" 


class Event {
public:
    int etype;
    double propensity;
     int site_id; 

    // 构造函数签名也必须更新
    // Event(int etype, double propensity, const Site& site); // 旧
    Event(int etype, double propensity, int site_id); 
};


#endif // EVENT_H