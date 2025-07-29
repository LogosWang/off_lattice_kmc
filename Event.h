#ifndef EVENT_H
#define EVENT_H

// 包含 Site 类的定义，因为事件的执行和倾向计算可能需要 Site 对象
#include "Site.h" 


class Event {
public:
    int etype;
    double propensity;
    const Site& site;
    Event(int etype,double propensity,const Site& site);
};


#endif // EVENT_H