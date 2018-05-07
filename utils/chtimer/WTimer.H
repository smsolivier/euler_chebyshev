#ifndef __WTIMER_H__
#define __WTIMER_H__
#include <chrono> 

using namespace std; 
class WTimer {
public:
	static WTimer& instance() {
		static WTimer wt; 
		return wt; 
	}
	~WTimer(); 
	WTimer(); 
private:
	chrono::time_point<chrono::system_clock> m_start; 
	chrono::duration<double> m_el; 
}; 

#endif