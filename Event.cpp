#include "Event.h" // Always include the corresponding header file
#include <cmath> 
Event::Event(int i, double prop, const Site& s) 
: etype(i), propensity(prop), site(s) {}