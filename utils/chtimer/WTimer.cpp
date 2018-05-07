#include "WTimer.H"
#include <iostream> 

WTimer::WTimer() {
	m_start = chrono::system_clock::now(); 
}

WTimer::~WTimer() {
	m_el = chrono::system_clock::now() - m_start; 
	cout << "Wall Time = " << m_el.count() << " seconds" << endl; 
}