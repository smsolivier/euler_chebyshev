#include "WTimer.H"
#include <iostream> 

WTimer::WTimer() {
	m_start = chrono::system_clock::now(); 
}

WTimer::~WTimer() {
	if (getenv("NWTIMER") != NULL) return; 
	m_el = chrono::system_clock::now() - m_start; 
	if (m_el.count() > 60) {
		cout << "Wall Time = " << m_el.count()/60 << " minutes" << endl; 
	} else {
		cout << "Wall Time = " << m_el.count() << " seconds" << endl; 
	}
}