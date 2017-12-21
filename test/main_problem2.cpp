# include <iostream>
# include "../rhsOdeProblem.H" 
# include <iomanip>
# include <cmath>
# include "../Euler/ForwardEulerSolver.H"
# include "../Euler/BackwardEulerSolver.H"
# include "../RungeKutta/ModifiedEuler/ModifiedEulerSolver.H"
# include "../RungeKutta/Heun/HeunSolver.H"
# include "../RungeKutta/RungeKutta4th/RungeKutta4Solver.H"
# include "../MultiStep/LeapFrogSolver.H"
# include "../RungeKutta/CrankNicholson/CrankNicholsonSolver.H"
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
using namespace mg::numeric::ode;
 
auto numFun = [](double t , double u) {return -20*u+20*sin(t)+cos(t) ; } ;

auto exacFun = [](double t, double y) { return exp(-20*t)+sin(t) ; };


int main(){
  
   const double t0 = 0.0;
   const double tf = 2.5;
   const double dt = 0.005;
   const double u0 = 1.0;
      
   string fname = "analitical_2.out" ; 
   
   rhsOdeProblem<double> p1(numFun, exacFun , t0, tf , dt, u0 , fname);
         
   ForwardEulerSolver<double> feuler1(p1) ;
   feuler1.solve("fwdEuler_2.out") ;
   
   BackwardEulerSolver<double> beuler1(p1) ;
   beuler1.solve("bwdEuler_2.out") ;
/*
   ModifiedEulerSolver<double> meuler1(p1) ;
   meuler1.solve("ModEuler_2.out") ;
*/
   
   CrankNicholsonSolver<double> cn(p1);
   cn.solve("CrankNicholson_2.out");

      

   HeunSolver<double> heun1(p1);   
   heun1.solve("Heun_2.out"); 

   RungeKutta4Solver<double> rk4_1(p1);
   rk4_1.solve("RK4_2.out");   
  
   LeapFrogSolver<double> leapFrog(p1);
   leapFrog.solve("LeapFrog_2.out");
  
  return 0;    
}

