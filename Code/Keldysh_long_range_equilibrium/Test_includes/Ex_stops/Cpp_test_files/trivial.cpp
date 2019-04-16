#include <iostream> 
#include <stdio.h>
#include <stdio.h>
#include <string.h>

using namespace std;

class A{
	public: 
		int a;
		A(int a_in): a(a_in){}
		virtual int a_times_b(int b)=0;
};

class B: public A{
	public:
		using A::a;
		using A::A;
		int a_times_b(int b){return a*b;}
};

int main(){
	const int mode=2;
	int i=5;
	B b(i);
	cout<<b.a_times_b(2)<<endl;
	return 0;
}
