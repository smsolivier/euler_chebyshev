#pragma once

#include <complex> 
#include <iostream> 
#include "WTimer.H"
#include "CH_Timer.H"

#define ERROR(message) {\
	std::cerr << "ERROR in " << __func__ << " (" << __FILE__ \
		<< " " << __LINE__  << "): " << message << std::endl; \
	exit(0); }

#ifdef WARNINGS 
#define WARNING(message) {\
	std::cerr << "WARNING (" << __func__ << " in " << __FILE__ \
		<< " " << __LINE__ << "): " << message << std::endl; }
#else 
#define WARNING(message) 
#endif 

#if defined ErrorCheck && ErrorCheck > 1
#define CHECK(test, message) \
	if (!(bool)(test)) {\
	std::cerr << "CHECK FAILED in " << __func__ << " (" << __FILE__ \
		<< " " << __LINE__  << "): " << message << std::endl; \
	exit(0); \
	}
#elif defined ErrorCheck
#define CHECK(test, message) \
	if (!(bool)(test)) {\
	std::cerr << "CHECK FAILED in " << __func__ << " (" << __FILE__ \
		<< " " << __LINE__  << "): " << message << std::endl; \
	}
#else 
#define CHECK(bool, message) 
#endif

#ifdef WARNINGS 
#define WCHECK(test, message) \
	if (!(bool)(test)) {\
	std::cerr << "WARNING: " << message << " (" << __func__ << " in " << __FILE__ \
		<< " " << __LINE__  << "): " << std::endl; \
	}
#else 
#define WCHECK(test, message) 
#endif

typedef std::complex<double> cdouble; 

// start wall timer 
static WTimer& wt = WTimer::instance();