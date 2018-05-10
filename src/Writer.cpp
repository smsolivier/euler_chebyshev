#include "Writer.H"
#include "CH_Timer.H"
#include "VisitWriter.H"
#include "Common.H"

Writer::Writer(string name) {
#ifndef NWRITE
	m_name = name; 
	m_count = 0; 
	m_writes = 0; 
	m_f = 1; 
	m_out.open(name+".visit"); 
	m_out << "!NBLOCKS " << 1 << endl; 
#endif
}

Writer::~Writer() {
	m_out.close(); 
}

void Writer::add(Scalar& a_scalar, string a_name) {
	m_scalars.push_back(&a_scalar); 
	m_scalar_names.push_back(a_name); 
}

void Writer::add(Vector& a_vector, string a_name) {
	m_vectors.push_back(&a_vector); 
	m_vector_names.push_back(a_name); 
}

void Writer::setFreq(int a_f) {m_f = a_f; } 

void Writer::write(bool force) {
#ifndef NWRITE 
	CH_TIMERS("write to vtk"); 

	if (m_count++%m_f != 0 && !force) {
		return; 
	}

	int nvars = m_scalars.size() + m_vectors.size(); 
	int n_scalars = m_scalars.size(); 
	int n_vectors = m_vectors.size(); 

	array<int,DIM> dims; 
	if (n_scalars > 0) {
		dims = m_scalars[0]->dims(); 
	} else if (n_vectors > 0) {
		dims = m_vectors[0]->dims(); 
	} else {
		ERROR("variables list is empty"); 
	}

	// copy if not in physical space 
	vector<Scalar> scalars(n_scalars); 
	vector<Vector> vectors(n_vectors); 
	// perform inverse if necessary 
	for (int i=0; i<n_scalars; i++) {
		if (m_scalars[i]->isFFC()) {
			m_scalars[i]->inverse(scalars[i]); // out of place transform 
		} else {
			scalars[i] = *m_scalars[i]; // copy pointer over 
		}
	}
	for (int i=0; i<n_vectors; i++) {
		if (m_vectors[i]->isFFC()) {
			m_vectors[i]->inverse(vectors[i]); 
		} else {
			vectors[i] = *m_vectors[i]; 
		}
	}

	// deep copy to float** 
	float** vars = new float*[nvars]; 
	for (int i=0; i<n_scalars; i++) {
		vars[i] = new float[scalars[i].size()]; 
		for (int j=0; j<scalars[i].size(); j++) {
			vars[i][j] = scalars[i][j].real(); 
		}
	}

	// vector version 
	for (int i=0; i<n_vectors; i++) {
		int vind = i + n_scalars; 
		vars[vind] = new float[vectors[0][0].size()*DIM]; 
		for (int j=0; j<vectors[0][0].size(); j++) {
			for (int d=0; d<DIM; d++) {
				vars[vind][DIM*j+d] = vectors[i][d][j].real(); 
			}
		}
	}

	// number of dimensions for each variable 
	int vardim[nvars]; 
	for (int i=0; i<n_scalars; i++) {
		vardim[i] = 1; 
	}
	for (int i=n_scalars; i<nvars; i++) {
		vardim[i] = DIM; 
	}

	// centering of the variables 
	int centering[nvars]; 
	for (int i=00; i<nvars; i++) {
		centering[i] = 1; 
	}

	// variable names 
	const char* varnames[nvars]; 
	for (int i=0; i<n_scalars; i++) {
		varnames[i] = m_scalar_names[i].c_str(); 
	}
	for (int i=n_scalars; i<nvars; i++) {
		varnames[i] = m_vector_names[i-n_scalars].c_str(); 
	}

	// grid locations 
	vector<float> x(dims[0]); 
	vector<float> y(dims[1]); 
	vector<float> z(dims[2]); 

	double xb = 2*M_PI; 
	double yb = 2*M_PI; 
	for (int i=0; i<dims[0]; i++) {
		x[i] = i*xb/dims[0]; 
	}
	for (int i=0; i<dims[1]; i++) {
		y[i] = i*yb/dims[1]; 
	}
	for (int i=0; i<dims[2]; i++) {
		z[i] = cos(i*M_PI/(dims[2]-1)); 
	}

	// append to master file 
	m_out << m_name << m_writes << ".vtk" << endl; 
	string fname = m_name + to_string(m_writes++); 

	write_rectilinear_mesh(fname.c_str(), VISIT_ASCII, &dims[0], 
		&z[0], &y[0], &x[0], nvars, &vardim[0], &centering[0], varnames, vars);

	// clean up pointers 
	for (int i=0; i<nvars; i++) {
		delete[] vars[i]; 
	}

	delete vars; 
#endif
}