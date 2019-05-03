#include <iostream>
#include <string.h> 
#include <time.h> 

using namespace std;

template<typename T> class A{
	public:
		T a;
		A(T a_in):a(a_in){}
		T operator()(){ return a;}
};

template<typename U> class B{
	public:
		U &Ar;
		B(U &Ar_in):Ar(Ar_in){}
		typename std::result_of<U()>::type operator()(){return Ar();}
};

int main(){
	A<int> a(5);
	cout<<"a()="<<a()<<endl;
	B<A<int> > b(a);
	cout<<"b()="<<b()<<endl;

	return 0;
}
