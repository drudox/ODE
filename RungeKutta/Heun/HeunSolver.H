# ifndef __HEUN_SOLVER_H__
# define __HEUN_SOLVER_H__

# include "../RungeKutta.H"
# include "../../rhsOdeProblem.H"

namespace mg {
                namespace numeric {
                                     namespace ode {


/*-------------------------------------------------------------------|  
 *                                                                      
 *   @brief Solve the ODE problem :  y' = f(t,y)                     
 *                                                                   
 *   Runge-Kutta Scheme (2 order accuracy) > Heun Method            
 *                                                                   
 *   @author Marco Ghiani                                               
 *   @date Dec 2017                                                     
 *   @place Glasgow UK                                               
 *                                                                   
 -------------------------------------------------------------------*/



template <typename Type>
class HeunSolver 
                  : public RungeKutta<Type>
{
  

  public:  

      HeunSolver(const rhsOdeProblem<Type> & that) noexcept  :
                                                                RungeKutta<Type>{that} 
                  {}
      
      virtual ~HeunSolver() = default ;
      
      using OdeSolver<Type>::rhs;
      
      void solve(const std::string filename) override final  ;
      void solve() noexcept override final                   ;
     
    
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



template <typename Type>
inline void HeunSolver<Type>::solve(const std::string filename) {

      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {     
         std::cerr << "Error opening file " << filename << std::endl ;  
      }
      else
      {
         std::cout << "Running Heun (RK -2nd ord) Solver" << std::endl;
         
         t.at(0) = t0();
         u.at(0) = u0();
         
         f << t.at(0) << ' ' << u.at(0) << std::endl; 
         
         for(auto i=1; i <= Ns ; i++ )
         {
            t.at(i) = t.at(i-1) + dt() ;
            
            k1 = rhs.f(t.at(i-1),u.at(i-1));   
            k2 = rhs.f(t.at(i-1)+dt() , u.at(i-1)+dt()*k1); 
            
            u.at(i) = u.at(i-1) + dt()/2 *(k1+k2) ;  
            f << t.at(i) << ' ' << u.at(i) << std::endl ;
         } 
         std::cout << "... Done " << std::endl;  
      }
      f.close();
}

template <typename Type>
inline void HeunSolver<Type>::solve() noexcept
{
     std::cout << "Running Heun (RK -2nd ord) Solver" << std::endl;
       

       t.at(0) = t0();
       u.at(0) = u0();
         
       std::cout << t.at(0) << ' ' << u.at(0) << std::endl; 
         
      for(auto i=1; i <= Ns ; i++ )
      {
          t.at(i) = t.at(i-1) + dt() ;
          k1 = rhs.f(t.at(i-1),u.at(i-1));   
          k2 = rhs.f(t.at(i-1)+dt() , u.at(i-1)+dt()*k1); 
            
          u.at(i) = u.at(i-1) + dt()/2 *(k1+k2) ;  
          std::cout << t.at(i) << ' ' << u.at(i) << std::endl ;
      } 
      std::cout << "... Done " << std::endl;  
}
  
  }//ode
 }//numeric
}//mg 
# endif
