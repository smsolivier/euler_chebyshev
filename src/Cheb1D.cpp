#include "Cheb1D.H"

Cheb1D::Cheb1D() {
	m_plan = NULL; 
}

Cheb1D::~Cheb1D() {
	if (m_plan != NULL) fftw_destroy_plan(m_plan); 
}

Cheb1D::Cheb1D(int N, int stride, cdouble* input) {
	init(N, stride, input); 
}

void Cheb1D::init(int N, int stride, cdouble* input) {
	double* in = reinterpret_cast<double*>(input); 
	m_N = N; 
	m_stride = stride; 
	m_stride *= 2; // double to account for complex part 
#ifdef MEASURE 
	int measure = FFTW_MEASURE; 
#else 
	int MEASURE = FFTW_ESTIMATE; 
#endif
	fftw_r2r_kind type = FFTW_REDFT00; // DCT type 1
	// setup plan
	m_plan = fftw_plan_many_r2r(
		1, // rank of transform 
		&m_N, // size of transform 
		1, // number of transforms 
		in, // input array pointer 
		NULL, // inembed 
		m_stride, // input stride 
		0, // idist 
		in, // output array pointer 
		NULL, // onembed 
		m_stride, // output stride 
		0, // odist 
		&type, // type, DCT 1 
		measure // FFTW flags, measure or estimate 
		); 
}

void Cheb1D::transform(cdouble* input, int DIR) {

	double* cast = reinterpret_cast<double*>(input); 

	// forward transform to chebyshev space 
	if (DIR == -1) {
		fftw_execute_r2r(m_plan, cast, cast); 
		for (int i=0; i<m_N; i++) {
			cast[i*m_stride] /= m_N-1; 
		}

		cast[0] /= 2; 
	}

	// inverse transform to physical space 
	else if (DIR == 1) {
		// denormalize 
		for (int i=0; i<m_N; i++) {
			cast[i*m_stride] *= m_N-1; 
		}

		input[0] *= 2; 
		fftw_execute_r2r(m_plan, cast, cast); 
		for (int i=0; i<m_N; i++) {
			cast[i*m_stride] /= 2*(m_N-1); 
		}
	}

	else {
		ERROR("DIR must be -1 or 1"); 
	}
}