# include <iostream>
# include <iomanip>
# include <string>
# include "../rhsOdeProblem.H" 
# include <cmath>
# include "../Euler/ForwardEulerSolver.H"
# include "../Euler/BackwardEulerSolver.H"
# include "../RungeKutta/ModifiedEuler/ModifiedEulerSolver.H"
# include "../RungeKutta/Heun/HeunSolver.H"
# include "../RungeKutta/RungeKutta4th/RungeKutta4Solver.H"
# include "../MultiStep/LeapFrogSolver.H"
# include "../MultiStep/AdamsMethods/AdamsBashforth/AdamsBashforth2ndSolver.H"
# include "../MultiStep/AdamsMethods/AdamsBashforth/AdamsBashforth3thSolver.H"
# include "../MultiStep/AdamsMethods/AdamsBashforth/AdamsBashforth4thSolver.H"
# include "../MultiStep/AdamsMethods/AdamsBashforth/AdamsBashforth5thSolver.H"
# include "../MultiStep/AdamsMethods/AdamsMoulton/AdamsMoulton2ndSolver.H"
# include "../MultiStep/AdamsMethods/AdamsMoulton/AdamsMoulton3thSolver.H"
# include "../MultiStep/AdamsMethods/AdamsMoulton/AdamsMoulton4thSolver.H"
# include "../MultiStep/AdamsMethods/AdamsMoulton/AdamsMoulton5thSolver.H"



using namespace std;
using namespace mg::numeric::ode ;
 
//auto numFun(double t, double u) { return -10*(t-1)*u; }
auto numFun = [](double t, double u) { return -2.0*t*u*u; } ;  //f(t,u)

auto exacFun = [](double t, double y) { return 1.0/(1+pow(t,2)) ; } ;   // y(t)


int main(){
  
   const double t0 = -5.0;
   const double tf = 5.0;
   const double dt = 0.04;
   const double u0 = 1.0/26.0 ;  //exp(-5.0);
      
   string fname = "analitical_4.out" ; 
   
   rhsOdeProblem<double> p1(numFun, exacFun , t0, tf , dt, u0 , fname);
         
   ForwardEulerSolver<double> feuler1(p1) ;
   feuler1.solve("fwdEuler_4.out") ;
   
   BackwardEulerSolver<double> beuler1(p1) ;
   beuler1.solve("bwdEuler_4.out") ;

   ModifiedEulerSolver<double> meuler1(p1) ;
   meuler1.solve("ModEuler_4.out") ;

   HeunSolver<double> heun1(p1);   
   heun1.solve("Heun_4.out"); 

   RungeKutta4Solver<double> rk4_1(p1);
   rk4_1.solve("RK4_4.out");   

   LeapFrogSolver<double> leapFrog(p1);
   leapFrog.solve("LeapFrog_4.out");
   
   AdamsBashforth2ndSolver<double> ab2(p1);
   ab2.solve("AdamBashforth2_4.out");
 
   AdamsBashforth3thSolver<double> ab3(p1);
   ab3.solve("AdamBashforth3_4.out");

   AdamsBashforth4thSolver<double> ab4(p1);
   ab4.solve("AdamBashforth4_4.out");

   AdamsBashforth5thSolver<double> ab5(p1);
   ab5.solve("AdamBashforth5_4.out");
 
   AdamsMoulton2ndSolver<double> am2(p1);
   am2.solve("AdamMoulton2_4.out");
 
   AdamsMoulton3thSolver<double> am3(p1);
   am3.solve("AdamMoulton3_4.out");

   AdamsMoulton4thSolver<double> am4(p1);
   am4.solve("AdamMoulton4_4.out");

   AdamsMoulton5thSolver<double> am5(p1);
   am5.solve("AdamMoulton5_4.out");



  return 0;    
}
