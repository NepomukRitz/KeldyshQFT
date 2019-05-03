/*
 * compute self-energy of the SIAM with a static vertex
 */

#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <array>

const std::complex<double> I (0,1);

// #include "matrix.h"
#include <Eigen/Dense> // probably not necessary
// #include <mach/mach.h> // memory usage, TODO: figure out how to check memory

#include "data_structures.h"   // defines essential data types
#include "util.h"              // time
#include "physics_functions.h" // functions providing physical objects

using namespace std;
using Eigen::MatrixXcd;

/* // attempt to use the "matrix" class, maybe better use "vector"
class SelfEnergy{
  public:
    matrix<matrix<complex<double> > > > SE_up;
    matrix<matrix<complex<double> > > > SE_down;
}

SelfEnergy::SelfEnergy(int U, int Nfreq) {
  for (i=0, i<Nfreq, ++i) {
    SE_up(0,0)
  }

}
*/

/* // attempt to use the "eigen" library, maybe better use "vector"
class SelfEnergy0{
  public:
    MatrixXcd SE(2, 2, 2, 1);
    SelfEnergy0(int N_omega);
    void initialize(double U);
};

SelfEnergy0::SelfEnergy0(int N_omega) {
  SE.resize(2, 2, 2, N_omega);
}

void SelfEnergy0::initialize(double U) {
  int N_omega = SE(0, 0, 0).size();
  for (int i_sigma=0; i_sigma<2; ++i_sigma) {
    for (int i_omega=0; i_omega<N_omega; ++i_omega) {
      SE(i_sigma, 0, 1, i_omega) = U / 2;  // retarded self-energy
      SE(i_sigma, 1, 0, i_omega) = U / 2;  // advanced self-energy
    }
  }
}
*/




int main() {

  double t0 = get_time(); // set start time

  /* // MEMORY USAGE (??) // TODO: check this, especially for Linux
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

  if (KERN_SUCCESS != task_info(mach_task_self(),
                                TASK_BASIC_INFO, (task_info_t)&t_info,
                                &t_info_count))
  {
    return -1;
  }

  cout << t_info.resident_size << "    " << t_info.virtual_size << endl;
  */


  /* // PLAY AROUND
  cout << "Hello" << endl;

  comp test_c (1, 1);
  vector<comp> test_v (3) ;
  vector<vector<comp> > test_vv (3, vector<comp> (3));

  cout << test_c << endl;
  cout << test_v.size() << endl;
  cout << test_vv.size() << endl;
  cout << test_vv[0].size() << endl;

  test_c = comp(2, 2);
  // test_v2 = test_v[1];
  test_vv[0][0] = 1;

  cout << test_c << endl;
  cout << test_vv[0][0] << endl;

  vector<comp> test_v2;
  test_v2 = vector<comp> (2);
  test_v2[0] = 2;

  cout << endl << "......" << endl << endl;

  cvec test (3);
  comp a = test.front();
  cout << test.front() << endl;
  cout << test.back() << endl;

  */


  /* // PLAY AROUND WITH SELF-ENERGY

  SelfEnergy SE (10);
  SE.initialize(1);

  cout << SE.SE[1][0][1][2] << endl;

  int a0[] = {1, 2, 3};
  rvec a (a0, a0 + sizeof(a0)/ sizeof(a0[0]));
  double mu = 0.5;
  double T = 1.5;

  rvec b = fermiFunction(a, T, mu);

  cout << 1./(exp((a[0] - mu)/T) + 1.) << endl;

  for (comp x : b) {
    cout << x << ",";
  }
  cout << endl;

  */

  cout << "----------" << endl;

  /* // TRY OUT KELDYSH MATRIX CLASS WITH INTEGER ARGUMENT

  KeldyshM2a<comp> m1 (0);
  KeldyshM2a<comp> m2 (0);
  cout << m1(0, 0) << endl;
  m1(0, 0) = 1;
  m1(0, 1) = 1;

  m2(0, 1) = 1;
  m2(1, 1) = 1;
  cout << m1(0, 0) << m1(0, 1) << m1(1, 0) << m1(1, 1) << endl;
  cout << m2(0, 0) << m2(0, 1) << m2(1, 0) << m2(1, 1) << endl;

  KeldyshM2a<comp> m (0);
  m = m1 * m2;
  cout << m(0, 0) << m(0, 1) << m(1, 0) << m(1, 1) << endl;

  */

  /* // TRY OUT KELDYSH MATRIX CLASS WITH ARRAY/VECTOR ARGUMENT

  KeldyshM2a<carray<3> > m3; // ({0, 0, 0});
  KeldyshM2a<carray<3> > m4; // ({0, 0, 0});
  //KeldyshM2a<cvec> m3 ({0, 0, 0});
  //KeldyshM2a<cvec> m4 ({0, 0, 0});

  m3(0, 0) = {1, 2, 3};
  m3(0, 1) = {1, 2, 3};
  m4(0, 1) = {1, 2, 3};
  m4(1, 1) = {1, 2, 3};

  KeldyshM2a<carray<3> > m5; // ({0, 0, 0});
  //KeldyshM2a<cvec> m5 ({0, 0, 0});
  m5 = m3 * m4;

  cout << m5(0, 1)(2) << endl;

  KeldyshM2<cvec> m6 (3);
  KeldyshM2<cvec> m7 ({0, 0, 0});

  m6(0, 0) = {1, 2, 3};
  m6(0, 1) = {1, 2, 3};
  m7(0, 1) = {1, 2, 3};
  m7(1, 1) = {1, 2, 3};

  KeldyshM2<cvec> m8;

  m8 = m6 * m7;

  cout << m8(0, 1)(2) << endl;

  // */


  //vector<KeldyshM2<cvec> > SE (N_omega);
  //cout << SE[0](0, 0)(0) << endl;

  V2P SE1 (3);
  cout << SE1(0, 0, 0, 1) << endl;

  V2P SE2 (3);
  V2P SE3 (3);
  SE3 = SE1 * SE2;
  cout << SE3(0, 0, 0, 2) << endl;



  int N_omega = 11;
  double omega_1 = -1;
  double omega_2 = 1;

  rvec eps_sigma ({-1, 1});

  const double U = 1;
  const double T = 1;
  double mu = 0;
  double Gamma = 1;
  double Lambda = 1;

  rvec omega (N_omega);
  for (int i=0; i<N_omega; ++i) {
    omega(i) = omega_1 + i*(omega_2-omega_1)/(N_omega-1);
  }

  V2P G0 (N_omega);
  G0 = get_G0(omega, eps_sigma, T, mu, Gamma, Lambda);

  V2P SE = initialize_SE(N_omega, U);



  // usleep(1000); // wait for a certain number of microseconds
  get_time(t0); // print time elapsed since time t0

}
