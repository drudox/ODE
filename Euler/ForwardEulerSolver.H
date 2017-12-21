# ifndef __FORWARD_EULER_SOLVER_H__
# define __FORWARD_EULER_SOLVER_H__

# include "Euler.H"
# include "../rhsOdeProblem.H"

namespace mg {
                namespace numeric {
                                     namespace ode {


/*-----------------------------------------------------------------------------
 *    Compute Explicit Euler to solve a given (ODE) RHS problem 
 *    
 *
 *    @author Marco Ghiani 
 *    @place Dec 2017, Glasgow UK
 *
 -----------------------------------------------------------------------------*/




template<typename Type = double>
class ForwardEulerSolver :
                            public Euler<Type> 
{
      
    public:  
      ForwardEulerSolver(const rhsOdeProblem<Type> & that) noexcept :
                                                                        Euler<Type>{that} 
                  {}
      
      virtual ~ForwardEulerSolver() = default ;

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

};

//------------------  Implementation (to be put into .cpp file)   -----------------  //


template<typename Type>
inline void ForwardEulerSolver<Type>::solve(const std::string filename) {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {     
         throw std::runtime_error(">> Error opening File in FWD euler solver <<" );
      }
      else
      {
         std::cout << "Running ForwardEuler Solver" << std::endl;
      
         t.at(0) = t0();
         u.at(0) = u0();
         
         f << t.at(0) << ' ' << u.at(0) << std::endl; 
         
         for(auto i=1; i <= Ns ; i++ )
         {
            t.at(i) = t.at(i-1) + dt() ;
            u.at(i) = u.at(i-1) + dt() *rhs.f(t.at(i-1), u.at(i-1)) ;  
            f << t.at(i) << ' ' << u.at(i) << std::endl ;
         } 
         std::cout << "... Done " << std::endl;  
      
         f.close();
      } 
}


template<typename Type>
inline void ForwardEulerSolver<Type>::solve() noexcept 
{
     std::cout << "Running ForwardEuler Solver" << std::endl;
      
       t.at(0) = t0();
       u.at(0) = u0();
         
       std::cout << t.at(0) << ' ' << u.at(0) << std::endl; 
         
      for(auto i=1; i <= Ns ; i++ )
      {
          t.at(i) = t.at(i-1) + dt() ;
          u.at(i) = u.at(i-1) + dt() *rhs.f(t.at(i-1), u.at(i-1)) ;  
          std::cout << t.at(i) << ' ' << u.at(i) << std::endl ;
      } 
      std::cout << "... Done " << std::endl;  
}
  
  }//ode
 }//numeric
}//mg
# endif
