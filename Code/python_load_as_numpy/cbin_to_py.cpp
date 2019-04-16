#include <boost/python.hpp>
#include <iostream>
#include "matrix.h"
#include <iomanip>
#include <ctime>
#include <math.h>
#include <string.h>
#include <complex> 

using namespace boost::python;

//First level structures:

struct Mi{
	matrix<int> data;
 	int load_components(const char *folder ,char *v_name){ 
		data.load(folder,v_name);
		return 0;
	}
 	int save_components(const char *folder ,char *v_name){ 
		data.save(folder,v_name);
		return 0;
	}
	void resize(int N1, int N2){
		data.resize(N1,N2);
	}
 	void set_component(int i1, int i2, int x){ 
		data(i1,i2) = x;
	}
 	int component(int i1, int i2){ 
		return data(i1,i2);
	}

};

struct Md{
	matrix<double> data;
 	int load_components(const char *folder ,char *v_name){ 
		data.load(folder,v_name);
		return 0;
	}
 	int save_components(const char *folder ,char *v_name){ 
		data.save(folder,v_name);
		return 0;
	}
	void resize(int N1, int N2){
		data.resize(N1,N2);
	}
 	void set_component(int i1, int i2, double x){ 
		data(i1,i2) = x;
	}
	
 	double component(int i1, int i2){ 
		return data(i1,i2);
	}

};

struct Mcd{
	matrix<std::complex<double> > data;
 	int load_components(const char *folder ,char *v_name){ 
		data.load(folder,v_name);
		return 0;
	}
 	int save_components(const char *folder ,char *v_name){ 
		data.save(folder,v_name);
		return 0;
	}
	void resize(int N1, int N2){
		data.resize(N1,N2);
	}
 	void set_component(int i1, int i2, std::complex<double> x){ 
		data(i1,i2) = x;
	}
 	std::complex<double> component(int i1, int i2){ 
		return data(i1,i2);
	}

};

struct Si{
	syma<int> data;
 	int load_components(const char *folder ,char *v_name){ 
		data.load(folder,v_name);
		return 0;
	}
 	int save_components(const char *folder ,char *v_name){ 
		data.save(folder,v_name);
		return 0;
	}
	void resize(int N1){
		data.resize(N1);
	}
 	void set_component(int i1, int i2, int x){ 
		if(i1 >=i2){
			data(i1,i2) = x;
		}
		else{
			data(i2,i1) =x;
		}
	}
 	int component(int i1, int i2){ 
	 	if(i1 >= i2){
		 	return data(i1,i2);
		}
		else{
		 	return data(i2,i1);
		}
	}

};

struct Sd{
	syma<double> data;
 	int load_components(const char *folder ,char *v_name){ 
		data.load(folder,v_name);
		return 0;
	}
 	int save_components(const char *folder ,char *v_name){ 
		data.save(folder,v_name);
		return 0;
	}
	void resize(int N1){
		data.resize(N1);
	}
 	void set_component(int i1, int i2, double x){ 
		if(i1 >=i2){
			data(i1,i2) = x;
		}
		else{
			data(i2,i1) =x;
		}
	}
 	double component(int i1, int i2){ 
	 	if(i1 >= i2){
		 	return data(i1,i2);
		}
		else{
		 	return data(i2,i1);
		}
	}

};

struct Scd{
	syma<std::complex<double> > data;
 	int load_components(const char *folder ,char *v_name){ 
		data.load(folder,v_name);
		return 0;
	}
 	int save_components(const char *folder ,char *v_name){ 
		data.save(folder,v_name);
		return 0;
	}
	void resize(int N1){
		data.resize(N1);
	}
 	void set_component(int i1, int i2, std::complex<double> x){ 
		if(i1 >=i2){
			data(i1,i2) = x;
		}
		else{
			data(i2,i1) =x;
		}
	}
	
 	std::complex<double> component(int i1, int i2){ 
	 	if(i1 >= i2){
		 	return data(i1,i2);
		}
		else{
		 	return data(i2,i1);
		}
	}

};

//Second level structures:

struct MMd{
	matrix<matrix<double> > data;
 	int load_components(const char *folder ,char *v_name){ 
		data.load(folder,v_name);
		return 0;
	}
 	int save_components(const char *folder ,char *v_name){ 
		data.save(folder,v_name);
		return 0;
	}
	void resize(int N1, int N2, int N3, int N4){
		data.resize(N1,N2);
		for(int i1=0; i1<N1; ++i1){
			for(int i2=0; i2<N2; ++i2){
				data(i1,i2).resize(N3,N4);
			}
		}
	}
 	void set_component(int i1, int i2, int i3, int i4, double x){ 
		data(i1,i2)(i3,i4) = x;
	}
	
 	double component(int i1, int i2, int i3, int i4){ 
		return data(i1,i2)(i3,i4);
	}

};


struct MMcd{
	matrix<matrix<std::complex<double> > > data;
 	int load_components(const char *folder ,char *v_name){ 
		data.load(folder,v_name);
		return 0;
	}
 	int save_components(const char *folder ,char *v_name){ 
		data.save(folder,v_name);
		return 0;
	}
	void resize(int N1, int N2, int N3, int N4){
		data.resize(N1,N2);
		for(int i1=0; i1<N1; ++i1){
			for(int i2=0; i2<N2; ++i2){
				data(i1,i2).resize(N3,N4);
			}
		}
	}
 	void set_component(int i1, int i2, int i3, int i4, std::complex<double> x){ 
		data(i1,i2)(i3,i4) = x;
	}
	
 	std::complex<double> component(int i1, int i2, int i3, int i4){ 
		return data(i1,i2)(i3,i4);
	}

};

struct MScd{
	matrix<syma<std::complex<double> > > data;
 	int load_components(const char *folder ,char *v_name){ 
		data.load(folder,v_name);
		return 0;
	}
 	int save_components(const char *folder ,char *v_name){ 
		data.save(folder,v_name);
		return 0;
	}
	void resize(int N1, int N2, int N3){
		data.resize(N1,N2);
		for(int i1=0; i1<N1; ++i1){
			for(int i2=0; i2<N2; ++i2){
				data(i1,i2).resize(N3);
			}
		}
	}
 	void set_component(int i1, int i2, int i3, int i4, std::complex<double> x){ 
		data(i1,i2)(i3,i4) = x;
	}
	
 	std::complex<double> component(int i1, int i2, int i3, int i4){ 
	 	if(i3 >= i4){
		 	return data(i1,i2)(i3,i4);
		}
		else{
		 	return data(i1,i2)(i4,i3);
		}
	}

};

//Third level structures:

struct MMMd{
	matrix<matrix<matrix<double> > > data;
 	int load_components(const char *folder ,char *v_name){ 
		data.load(folder,v_name);
		return 0;
	}
 	int save_components(const char *folder ,char *v_name){ 
		data.save(folder,v_name);
		return 0;
	}
	
	void resize(int N1, int N2, int N3, int N4, int N5, int N6){
		data.resize(N1,N2);
		for(int i1=0; i1<N1; ++i1){
			for(int i2=0; i2<N2; ++i2){
				data(i1,i2).resize(N3,N4);
				for(int i3=0; i3<N3; ++i3){
					for(int i4=0; i4<N4; ++i4){
						data(i1,i2)(i3,i4).resize(N5,N6);
					}
				}
			}
		}
	}
 	void set_component(int i1, int i2, int i3, int i4, int i5, int i6, double x){ 
		data(i1,i2)(i3,i4)(i5,i6) = x;
	}
	
 	double component(int i1, int i2, int i3, int i4, int i5, int i6){ 
		return data(i1,i2)(i3,i4)(i5,i6);
	}

};

struct MMMcd{
	matrix<matrix<matrix<std::complex<double> > > > data;
 	int load_components(const char *folder ,char *v_name){ 
		data.load(folder,v_name);
		return 0;
	}
 	int save_components(const char *folder ,char *v_name){ 
		data.save(folder,v_name);
		return 0;
	}
	
	void resize(int N1, int N2, int N3, int N4, int N5, int N6){
		data.resize(N1,N2);
		for(int i1=0; i1<N1; ++i1){
			for(int i2=0; i2<N2; ++i2){
				data(i1,i2).resize(N3,N4);
				for(int i3=0; i3<N3; ++i3){
					for(int i4=0; i4<N4; ++i4){
						data(i1,i2)(i3,i4).resize(N5,N6);
					}
				}
			}
		}
	}
 	void set_component(int i1, int i2, int i3, int i4, int i5, int i6, std::complex<double> x){ 
		data(i1,i2)(i3,i4)(i5,i6) = x;
	}
	
 	std::complex<double> component(int i1, int i2, int i3, int i4, int i5, int i6){ 
		return data(i1,i2)(i3,i4)(i5,i6);
	}

};
 	
struct MMScd{
	matrix<matrix<syma<std::complex<double> > > > data;
 	int load_components(const char *folder ,char *v_name){ 
		data.load(folder,v_name);
		return 0;
	}
 	int save_components(const char *folder ,char *v_name){ 
		data.save(folder,v_name);
		return 0;
	}
	void resize(int N1, int N2, int N3, int N4, int N5){
		data.resize(N1,N2);
		for(int i1=0; i1<N1; ++i1){
			for(int i2=0; i2<N2; ++i2){
				data(i1,i2).resize(N3,N4);
				for(int i3=0; i3<N3; ++i3){
					for(int i4=0; i4<N4; ++i4){
						data(i1,i2)(i3,i4).resize(N5);
					}
				}
			}
		}
	}
 	void set_component(int i1, int i2, int i3, int i4, int i5, int i6, std::complex<double> x){ 
		data(i1,i2)(i3,i4)(i5,i6) = x;
	}
	
 	std::complex<double> component(int i1, int i2, int i3, int i4, int i5, int i6){ 
	 	if(i5 >= i6){
		 	return data(i1,i2)(i3,i4)(i5,i6);
		}
		else{
		 	return data(i1,i2)(i3,i4)(i6,i5);
		}
	}

};

BOOST_PYTHON_MODULE(cbin_to_py)
{
//First level structures:

	class_<Mi>("Mi_py")
	      .def("load_components", &Mi::load_components)
	      .def("save_components", &Mi::save_components)
	      .def("resize", &Mi::resize)
	      .def("set_component", &Mi::set_component)
	      .def("component", &Mi::component)
	;
	
	class_<Md>("Md_py")
	      .def("load_components", &Md::load_components)
	      .def("save_components", &Md::save_components)
	      .def("resize", &Md::resize)
	      .def("set_component", &Md::set_component)
	      .def("component", &Md::component)
	;

	class_<Mcd>("Mcd_py")
	      .def("load_components", &Mcd::load_components)
	      .def("save_components", &Mcd::save_components)
	      .def("resize", &Mcd::resize)
	      .def("set_component", &Mcd::set_component)
	      .def("component", &Mcd::component)
	;

	class_<Si>("Si_py")
	      .def("load_components", &Si::load_components)
	      .def("save_components", &Si::save_components)
	      .def("resize", &Si::resize)
	      .def("set_component", &Si::set_component)
	      .def("component", &Si::component)
	;
	
	class_<Sd>("Sd_py")
	      .def("load_components", &Sd::load_components)
	      .def("save_components", &Sd::save_components)
	      .def("resize", &Sd::resize)
	      .def("set_component", &Sd::set_component)
	      .def("component", &Sd::component)
	;

	class_<Scd>("Scd_py")
	      .def("load_components", &Scd::load_components)
	      .def("save_components", &Scd::save_components)
	      .def("resize", &Scd::resize)
	      .def("set_component", &Scd::set_component)
	      .def("component", &Scd::component)
	;

//Second level structures:

	class_<MMcd>("MMcd_py")
	      .def("load_components", &MMcd::load_components)
	      .def("save_components", &MMcd::save_components)
	      .def("resize", &MMcd::resize)
	      .def("set_component", &MMcd::set_component)
	      .def("component", &MMcd::component)
	;
	
	class_<MMd>("MMd_py")
	      .def("load_components", &MMd::load_components)
	      .def("save_components", &MMd::save_components)
	      .def("resize", &MMd::resize)
	      .def("set_component", &MMd::set_component)
	      .def("component", &MMd::component)
	;
	
	class_<MScd>("MScd_py")
	      .def("load_components", &MScd::load_components)
	      .def("save_components", &MScd::save_components)
	      .def("resize", &MScd::resize)
	      .def("set_component", &MScd::set_component)
	      .def("component", &MScd::component)
	;

//Third level structures:
	class_<MMMd>("MMMd_py")
	      .def("load_components", &MMMd::load_components)
	      .def("save_components", &MMMd::save_components)
	      .def("resize", &MMMd::resize)
	      .def("set_component", &MMMd::set_component)
	      .def("component", &MMMd::component)
	;

	class_<MMMcd>("MMMcd_py")
	      .def("load_components", &MMMcd::load_components)
	      .def("save_components", &MMMcd::save_components)
	      .def("resize", &MMMcd::resize)
	      .def("set_component", &MMMcd::set_component)
	      .def("component", &MMMcd::component)
	;
	
	class_<MMScd>("MMScd_py")
	      .def("load_components", &MMScd::load_components)
	      .def("save_components", &MMScd::save_components)
	      .def("resize", &MMScd::resize)
	      .def("set_component", &MMScd::set_component)
	      .def("component", &MMScd::component)
	;

	
}


