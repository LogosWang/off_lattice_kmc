// Event.cpp

#include "Event.h" 
#include <cmath> 

// *** 关键修改：更新构造函数定义 ***
Event::Event(int i, double prop, int id) 
: etype(i), propensity(prop), site_id(id) {}

// 旧的构造函数定义 (应删除或注释掉):
// Event::Event(int i, double prop, const Site& s) 
// : etype(i), propensity(prop), site(s) {}