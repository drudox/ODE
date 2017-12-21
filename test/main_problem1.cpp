# include <iostream>
# include <iomanip>
# include <string>
# include "../rhsOdeProblem.H" 
# include <cmath>
# include "../Euler/ForwardEulerSolver.H"
# include "../Euler/BackwardEulerSolver.H"
//# include "ModifiedEulerSolver.H"
//# include "HeunSolver.H"
//# include "RungeKutta4Solver.H"
# include "../MultiStep/LeapFrogSolver.H"
//# include "ImplicitEulerSolver.H"
# include "../MultiStep/AdamsMethods/AdamsBashforth/AdamsBashforth2ndSolver.H"
# include "../MultiStep/AdamsMethods/AdamsBashforth/AdamsBashforth3thSolver.H"
# include "../MultiStep/AdamsMethods/AdamsBashforth/AdamsBashforth4thSolver.H"
# include "../MultiStep/AdamsMethods/AdamsBashforth/AdamsBashforth5thSolver.H"
# include "../MultiStep/AdamsMethods/AdamsMoulton/AdamsMoulton2ndSolver.H"
# include "../MultiStep/AdamsMethods/AdamsMoulton/AdamsMoulton3thSolver.H"
# include "../MultiStep/AdamsMethods/AdamsMoulton/AdamsMoulton4thSolver.H"
# include "../MultiStep/AdamsMethods/AdamsMoulton/AdamsMoulton5thSolver.H"
# include "../RungeKutta/CrankNicholson/CrankNicholsonSolver.H"

using namespace std;
using namespace mg::numeric::ode ;
 
//auto numFun(double t, double u) { return -10*(t-1)*u; }
auto numFun = [](double t, double u) { return -10*(t-1)*u; } ;

auto exacFun = [](double t, double y) { return exp(-5*pow((t-1),2) ); };


int main(){
  
   const double t0 = 0.0;
   const double tf = 2.0;
   const double dt = 0.025;
   const double u0 = exp(-5.0);
      
   string fname = "analitical_1.out" ; 
   
   rhsOdeProblem<double> p1(numFun, exacFun , t0, tf , dt, u0 , fname);
         
   ForwardEulerSolver<double> feuler1(p1) ;
   feuler1.solve("fwdEuler_1.out") ;
   
   BackwardEulerSolver<double> beuler1(p1) ;
   beuler1.solve("bwdEuler_1.out") ;
 
   
   AdamsBashforth2ndSolver<double> ab2(p1) ;
   ab2.solve("AdamBashforth2nd_1.out");   
   
   AdamsBashforth3thSolver<double> ab3(p1) ;
   ab3.solve("AdamBashforth3th_1.out");   

   AdamsBashforth4thSolver<double> ab4(p1) ;
   ab4.solve("AdamBashforth4th_1.out");   
 
   AdamsBashforth5thSolver<double> ab5(p1) ;
   ab5.solve("AdamBashforth5th_1.out");   

   AdamsMoulton2ndSolver<double> am2(p1) ;
   am2.solve("AdamMoulton2nd_1.out");   
   
   AdamsMoulton3thSolver<double> am3(p1) ;
   am3.solve("AdamMoulton3th_1.out");   
    
   AdamsMoulton4thSolver<double> am4(p1) ;
   am4.solve("AdamMoulton4th_1.out");   
   
   AdamsMoulton5thSolver<double> am5(p1) ;
   am5.solve("AdamMoulton5th_1.out");   
   
   CrankNicholsonSolver<double> cn(p1);
   cn.solve("CrankNicholson_1.out");

/*
   ModifiedEulerSolver<double> meuler1(p1) ;
   meuler1.solve("ModEuler_1.out") ;

   HeunSolver<double> heun1(p1);   
   heun1.solve("Heun_1.out"); 

   RungeKutta4Solver<double> rk4_1(p1);
   rk4_1.solve("RK4_1.out");   
*/
   LeapFrogSolver<double> leapFrog(p1);
   leapFrog.solve("LeapFrog_1.out");
   
  return 0;    
}
