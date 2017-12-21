# ifndef __RUNGEKUTTA4_SOLVER_H__
# define __RUNGEKUTTA4_SOLVER_H__

//# include "OdeSolver.H"
# include "../../rhsOdeProblem.H"
# include "../RungeKutta.H" 

namespace mg { 
                namespace numeric {
                                    namespace ode {


/*-------------------------------------------------------------------------------
 *    
 *    Perform Runge-Kutta (4th order accuracy) solution of (ODE) RHS problem
 *    dy/dt = f(y,t)
 *
 *    @Marco Ghiani October 2017, Glasgow UK
 *
 ------------------------------------------------------------------------------*/



template<typename Type= double>
class RungeKutta4Solver 
                         :   public  RungeKutta<Type>    
{
      
    public:  
      RungeKutta4Solver(const rhsOdeProblem<Type> & that) noexcept :
                                                                        RungeKutta<Type>{that} 
                  {}
      
      virtual ~RungeKutta4Solver() = default;
      

      using OdeSolver<Type>::rhs;

      void solve(const std::string filename) override final ;
      void solve() noexcept override final                     ;
     
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
      using RungeKutta<Type>::k3;
      using RungeKutta<Type>::k4;
 

};


template<typename Type>
inline void RungeKutta4Solver<Type>::solve(const std::string filename)  {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {     
         std::string mess = "Error opening file " + filename + " in RungeKutta4-Solver " ;  
         throw std::runtime_error(mess.c_str());   
      }
      else
      {
         std::cout << "Running Runge-Kutta 4th order Solver" << std::endl;
         

         t.at(0) = t0();
         u.at(0) = u0();
         
         f << t.at(0) << ' ' << u.at(0) << std::endl; 
         
         for(auto i=1; i <= Ns ; i++ )
         {
            t.at(i) = t.at(i-1) + dt() ;
            
            k1 = rhs.f(t.at(i-1)        , u.at(i-1)            );   
            k2 = rhs.f(t.at(i-1)+dt()/2 , u.at(i-1) + dt()/2*k1); 
            k3 = rhs.f(t.at(i-1)+dt()/2 , u.at(i-1) + dt()/2*k2);
            k4 = rhs.f(t.at(i-1)+dt()   , u.at(i-1) + dt()*k3  ); 

            u.at(i) = u.at(i-1) + dt()/6 *(k1+ 2*k2 + 2*k3 +k4) ;  
            f << t.at(i) << ' ' << u.at(i) << std::endl ;
         } 
         std::cout << "... Done " << std::endl;  
      
      f.close();
      }
}


template<typename Type>
inline void RungeKutta4Solver<Type>::solve() noexcept 
{
     std::cout << "Running Runge-Kutta 4th order Solver" << std::endl;
       
       t.at(0) = t0();
       u.at(0) = u0();
         
       std::cout << t.at(0) << ' ' << u.at(0) << std::endl; 
         
      for(auto i=1; i <= Ns ; i++ )
      {
          t.at(i) = t.at(i-1) + dt() ;

          k1 = rhs.f(t.at(i-1)        , u.at(i-1)            );   
          k2 = rhs.f(t.at(i-1)+dt()/2 , u.at(i-1) + dt()/2*k1); 
          k3 = rhs.f(t.at(i-1)+dt()/2 , u.at(i-1) + dt()/2*k2);
          k4 = rhs.f(t.at(i-1)+dt()   , u.at(i-1) + dt()*k3  ); 

          u.at(i) = u.at(i-1) + dt()/6 *(k1+ 2*k2 + 2*k3 +k4) ; 
          std::cout << t.at(i) << ' ' << u.at(i) << std::endl ;
      } 
      std::cout << "... Done " << std::endl;  
}

  }//ode
 }//numeric
}//mg 
# endif
