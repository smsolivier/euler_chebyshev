#pragma once

#include <complex> 
#include <iostream> 
#include "WTimer.H"
#include "CH_Timer.H"

#define ERROR(message) {\
	std::cout << "ERROR in " << __func__ << " (" << __FILE__ \
		<< " " << __LINE__  << "): " << message << std::endl; \
	exit(0); }

#ifndef NCHECK 
#define CHECK(test, message) \
	if (!(bool)test) {\
	std::cout << "CHECK FAILED in " << __func__ << " (" << __FILE__ \
		<< " " << __LINE__  << "): " << message << std::endl; \
	exit(0); }
#else 
#define CHECK(bool, message) 
#endif

typedef std::complex<double> cdouble; 

// start wall timer 
static WTimer& wt = WTimer::instance();