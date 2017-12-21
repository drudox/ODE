#ifndef __MODIFIED_EULER_H__
#define __MODIFIED_EULER_H__


# include "../../rhsOdeProblem.H"
# include "../RungeKutta.H"

namespace mg { 
               namespace numeric {
                                    namespace ode {

/*-----------------------------------------------------------------------
 *    compute the solution of RHS (ODE) problem  y' = f(t,y) 
 *    using (Runge Kutta 2th order accuracy) 
 *    Modified - explicit Euler       
 *    
 *    @Marco Ghiani Nov. 2017 Glasgow 
 *
 ----------------------------------------------------------------------*/




template <typename Type>
class ModifiedEulerSolver 
                          : public RungeKutta<Type>
{


   public:   

    ModifiedEulerSolver(const rhsOdeProblem<Type>& that) noexcept :
                                                                    RungeKutta<Type>{that} 
                    {}
    
    virtual ~ModifiedEulerSolver() = default;
    

    using OdeSolver<Type>::rhs;


    void solve(const std::string filename ) override final;
    void solve() noexcept override final                  ;

//--
   private:
      
      using OdeSolver<Type>::t  ;
      using OdeSolver<Type>::u  ;
      using OdeSolver<Type>::dt ; 
      using OdeSolver<Type>::t0 ;
      using OdeSolver<Type>::tf ;
      using OdeSolver<Type>::u0 ;
      
      using OdeSolver<Type>::Ns ;
      
      using OdeSolver<Type>::toll ;

      using RungeKutta<Type>::k1;
      using RungeKutta<Type>::k2;
      
};


//----------------- Implementation (to be put into .cpp file) ----------

template <typename Type>
inline void ModifiedEulerSolver<Type>::solve(const std::string fname)  
{
   std::ofstream f(fname , std::ios::out);

   if(!f)
   {
      std::string mess = "Error opening file " + fname + "in Modified Euler Solver " ; 
      throw std::runtime_error(mess.c_str());
   }
   
   t.at(0) = t0();
   u.at(0) = u0();   
   
   f << t.at(0) << ' ' << u.at(0) << std::endl;
   std::cout << "Running Modified Euler (RK -2nd ord) Solver" << std::endl;
   for(auto i=1 ; i <= Ns ; i++ )
   {
      t.at(i) = t.at(i-1) + dt() ;
      
      k1 = rhs.f(t.at(i-1) , u.at(i-1) );
      k2 = rhs.f(t.at(i-1)+ dt()/2 , u.at(i-1) + k1*dt()/2 );
      
      u.at(i) = u.at(i-1) + dt() * k2;
      
      f << t.at(i) << ' ' << u.at(i) << std::endl;
   }
   f.close();
   std::cout << "... Done" << std::endl;
}


template <typename Type>
inline void ModifiedEulerSolver<Type>::solve() noexcept 
{

   t.at(0) = t0();
   u.at(0) = u0();   
   
   std::cout << t.at(0) << ' ' << u.at(0) << std::endl;
   std::cout << "Running Modified (RK -2nd ord) Euler Solver" << std::endl;
   for(auto i=1 ; i <= Ns ; i++ )
   {
      t.at(i) = t.at(i-1) + dt() ;
      
      k1 = rhs.f(t.at(i-1) , u.at(i-1) );
      k2 = rhs.f(t.at(i-1)+ dt()/2 , u.at(i-1) + k1*dt()/2 );
      
      u.at(i) = u.at(i-1) + dt()* k2;
      
      std::cout << t.at(i) << ' ' << u.at(i) << std::endl;
   }
   std::cout << "... Done" << std::endl;
}
  
  }//ode
 }//numeric
}//mg 
# endif
