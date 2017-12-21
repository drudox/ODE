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

using namespace std;
using namespace mg::numeric::ode ;
 
//auto numFun(double t, double u) { return -10*(t-1)*u; }
auto numFun = [](double t, double u) { return t*u ; } ;

auto exacFun = [](double t, double y) { return exp(pow(t,2)/2.0); };


int main(){
  
   const double t0 = -2.0;
   const double tf = 2.0;
   const double dt = 0.3;
   const double u0 = exp(2.0);
      
   string fname = "analitical_3.out" ; 
   
   rhsOdeProblem<double> p1(numFun, exacFun , t0, tf , dt, u0 , fname);
         
   ForwardEulerSolver<double> feuler1(p1) ;
   feuler1.solve("fwdEuler_3.out") ;
   
   BackwardEulerSolver<double> beuler1(p1) ;
   beuler1.solve("bwdEuler_3.out") ;

   ModifiedEulerSolver<double> meuler1(p1) ;
   meuler1.solve("ModEuler_3.out") ;

   HeunSolver<double> heun1(p1);   
   heun1.solve("Heun_3.out"); 

   RungeKutta4Solver<double> rk4_1(p1);
   rk4_1.solve("RK4_3.out");   

   LeapFrogSolver<double> leapFrog(p1);
   leapFrog.solve("LeapFrog_3.out");
   
  return 0;    
}
