# ifndef __ADAMS_BASHFORTH_5STEP_SOLVER_H__
# define __ADAMS_BASHFORTH_5STEP_SOLVER_H__

# include "../AdamsMethods.H"
# include "../../../rhsOdeProblem.H"

namespace mg {
                namespace numeric {
                                     namespace ode {


/**----------------------------------------------------------------------------------------------
 * 
 * @class AdamsBashforth5thSolver
 * @brief Perform Explicit multi-step 4th order accuracy (4 step)
 *  solution of a given (ODE).  
 * 
 *
 * @author Marco Ghiani 
 * @date Dec 2017, Glasgow UK
 *
 *
 ---------------------------------------------------------------------------------------------------*/




template<typename Type = double>
class AdamsBashforth5thSolver :
                                 public AdamsMethods<Type> 
{
      
    public:  
      AdamsBashforth5thSolver(const rhsOdeProblem<Type> & that) noexcept :
                                                                             AdamsMethods<Type>{that} 
                  {}
      
      virtual ~AdamsBashforth5thSolver() = default ;

      using OdeSolver<Type>::rhs;

      void solve(const std::string filename) override final;
      void solve() noexcept override final                 ;
//
//
  private:

      using OdeSolver<Type>::t  ;
      using OdeSolver<Type>::u  ;
      using OdeSolver<Type>::dt ; 
      using OdeSolver<Type>::t0 ;
      using OdeSolver<Type>::tf ;
      using OdeSolver<Type>::u0 ;
      
      using OdeSolver<Type>::Ns ;

      using AdamsMethods<Type>::k1 ;
      using AdamsMethods<Type>::k2 ;
      using AdamsMethods<Type>::k3 ;
      using AdamsMethods<Type>::k4 ;


};

//------------------  Implementation (to be put into .cpp file)   -----------------  //


template<typename Type>
inline void AdamsBashforth5thSolver<Type>::solve(const std::string filename) {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {  
         std::string mess = "Error opening file: " + filename + " in Solver ADAMS-BASHFORTH 5 step" ;   
         throw std::runtime_error(mess.c_str());
      }
      else
      {
         std::cout << "Running Adams-Bashforth (5th order) Solver" << std::endl;
      
         t.at(0) = t0();
         u.at(0) = u0();   // initial Value 
         
         f << t.at(0) << ' ' << u.at(0) << std::endl; 
            
         // compute first point  (start-up the solver)  
         k1 = rhs.f(t.at(0)           , u.at(0)           );
         k2 = rhs.f(t.at(0)+ dt()/2.0 , u.at(0) + k1*dt()/2.0 );
         k3 = rhs.f(t.at(0)+ dt()/2.0 , u.at(0) + k2*dt()/2.0 );
         k4 = rhs.f(t.at(0)+ dt()     , u.at(0) + k3*dt() );
         
         u.at(1) =  u.at(0) + dt()/6.0 *(k1 + 2.*k2 + 2.*k3 + k4);  // RungeKutta 4th order 
   
         //  first point 
         t.at(1) = t.at(0) + dt() ;
   
         f << t.at(1) << ' ' << u.at(1) << std::endl; // write to file first point 

         // compute second point (start-up the solver)
         k1 = rhs.f(t.at(1)          , u.at(1)              );
         k2 = rhs.f(t.at(1)+ dt()/2. , u.at(1) + k1*dt()/2. );
         k3 = rhs.f(t.at(1)+ dt()/2. , u.at(1) + k2*dt()/2. );
         k4 = rhs.f(t.at(1)+ dt()    , u.at(1) + k3*dt()    );
         
         u.at(2) =  u.at(1) + dt()/6.0 *(k1 + 2.*k2 + 2.*k3 +k4);  // RungeKutta 4th order 
   
         //  second point  
         t.at(2) = t.at(1) + dt() ;
   
         f << t.at(2) << ' ' << u.at(2) << std::endl; // write to file second point 

         // compute third point (start-up the solver)
         k1 = rhs.f(t.at(2)          , u.at(2)              );
         k2 = rhs.f(t.at(2)+ dt()/2. , u.at(2) + k1*dt()/2. );
         k3 = rhs.f(t.at(2)+ dt()/2. , u.at(2) + k2*dt()/2. );
         k4 = rhs.f(t.at(2)+ dt()    , u.at(2) + k3*dt()    );
         
         u.at(3) =  u.at(2) + dt()/6.0 *(k1 + 2.*k2 + 2.*k3 +k4);  // RungeKutta 4th order 
   
         //  third point  
         t.at(3) = t.at(2) + dt() ;
   
         f << t.at(3) << ' ' << u.at(3) << std::endl; // write to file second point 
 
         // compute fourth point (start-up the solver)
         k1 = rhs.f(t.at(3)          , u.at(3)              );
         k2 = rhs.f(t.at(3)+ dt()/2. , u.at(3) + k1*dt()/2. );
         k3 = rhs.f(t.at(3)+ dt()/2. , u.at(3) + k2*dt()/2. );
         k4 = rhs.f(t.at(3)+ dt()    , u.at(3) + k3*dt()    );
         
         u.at(4) =  u.at(3) + dt()/6.0 *(k1 + 2.*k2 + 2.*k3 +k4);  // RungeKutta 4th order 
   
         //  fourth point  
         t.at(4) = t.at(3) + dt() ;
   
         f << t.at(4) << ' ' << u.at(4) << std::endl; // write to file fourth point 



         for(auto i=4; i < Ns ; i++ )
         {
            t.at(i+1) = t.at(i) + dt() ;
            u.at(i+1) = u.at(i) + dt()      *( 1901.0/720.0 * rhs.f(t.at(i)  , u.at(i)  )
                                              -1387.0/360.0 * rhs.f(t.at(i-1), u.at(i-1))
                                              + 109.0/30.0  * rhs.f(t.at(i-2), u.at(i-2)) 
                                              - 637.0/360.0 * rhs.f(t.at(i-3), u.at(i-3))  
                                              + 251.0/720.0 * rhs.f(t.at(i-4), u.at(i-4)) ) ;

            f << t.at(i+1) << ' ' << u.at(i+1) << std::endl ;
         } 
         std::cout << "... Done " << std::endl;  
      
         f.close();
      } 
}


template<typename Type>
inline void AdamsBashforth5thSolver<Type>::solve() noexcept 
{

     std::cout << "Running Adams-Bashforth (5th order) Solver" << std::endl;
      
         t.at(0) = t0();
         u.at(0) = u0();   // initial Value 
         
         std::cout << t.at(0) << ' ' << u.at(0) << std::endl; 
            
         // compute first point  (start-up the solver)  
         k1 = rhs.f(t.at(0)           , u.at(0)           );
         k2 = rhs.f(t.at(0)+ dt()/2.0 , u.at(0) + k1*dt()/2.0 );
         k3 = rhs.f(t.at(0)+ dt()/2.0 , u.at(0) + k2*dt()/2.0 );
         k4 = rhs.f(t.at(0)+ dt()     , u.at(0) + k3*dt() );
         
         u.at(1) =  u.at(0) + dt()/6.0 *(k1 + 2.*k2 + 2.*k3 + k4);  // RungeKutta 4th order 
   
         //  first point 
         t.at(1) = t.at(0) + dt() ;
   
         std::cout << t.at(1) << ' ' << u.at(1) << std::endl; // write to file first point 

         // compute second point (start-up the solver)
         k1 = rhs.f(t.at(1)          , u.at(1)              );
         k2 = rhs.f(t.at(1)+ dt()/2. , u.at(1) + k1*dt()/2. );
         k3 = rhs.f(t.at(1)+ dt()/2. , u.at(1) + k2*dt()/2. );
         k4 = rhs.f(t.at(1)+ dt()    , u.at(1) + k3*dt()    );
         
         u.at(2) =  u.at(1) + dt()/6.0 *(k1 + 2.*k2 + 2.*k3 +k4);  // RungeKutta 4th order 
   
         //  second point  
         t.at(2) = t.at(1) + dt() ;
   
         std::cout << t.at(2) << ' ' << u.at(2) << std::endl; // write to file second point 

         // compute third point (start-up the solver)
         k1 = rhs.f(t.at(2)          , u.at(2)              );
         k2 = rhs.f(t.at(2)+ dt()/2. , u.at(2) + k1*dt()/2. );
         k3 = rhs.f(t.at(2)+ dt()/2. , u.at(2) + k2*dt()/2. );
         k4 = rhs.f(t.at(2)+ dt()    , u.at(2) + k3*dt()    );
         
         u.at(3) =  u.at(2) + dt()/6.0 *(k1 + 2.*k2 + 2.*k3 +k4);  // RungeKutta 4th order 
   
         //  third point  
         t.at(3) = t.at(2) + dt() ;
   
         std::cout << t.at(3) << ' ' << u.at(3) << std::endl; // write to file second point 
 
         // compute fourth point (start-up the solver)
         k1 = rhs.f(t.at(3)          , u.at(3)              );
         k2 = rhs.f(t.at(3)+ dt()/2. , u.at(3) + k1*dt()/2. );
         k3 = rhs.f(t.at(3)+ dt()/2. , u.at(3) + k2*dt()/2. );
         k4 = rhs.f(t.at(3)+ dt()    , u.at(3) + k3*dt()    );
         
         u.at(4) =  u.at(3) + dt()/6.0 *(k1 + 2.*k2 + 2.*k3 +k4);  // RungeKutta 4th order 
   
         //  fourth point  
         t.at(4) = t.at(3) + dt() ;
   
         std::cout << t.at(4) << ' ' << u.at(4) << std::endl; // write to file fourth point 



         for(auto i=4; i < Ns ; i++ )
         {
            t.at(i+1) = t.at(i) + dt() ;
            u.at(i+1) = u.at(i) + dt()      *( 1901.0/720.0 * rhs.f(t.at(i)  , u.at(i)  )
                                              -1387.0/360.0 * rhs.f(t.at(i-1), u.at(i-1))
                                              + 109.0/30.0  * rhs.f(t.at(i-2), u.at(i-2)) 
                                              - 637.0/360.0 * rhs.f(t.at(i-3), u.at(i-3))  
                                              + 251.0/720.0 * rhs.f(t.at(i-4), u.at(i-4)) ) ;

            std::cout << t.at(i+1) << ' ' << u.at(i+1) << std::endl ;
         } 
         std::cout << "... Done " << std::endl; 

 
}
  
  }//ode
 }//numeric
}//mg
# endif
