# ifndef __EULER_ABSTRACT_INTERFACE_H__
# define __EULER_ABSTRACT_INTERFACE_H__

# include "../OdeSolver.H"
# include "../rhsOdeProblem.H"

namespace mg {
                namespace numeric {
                                     namespace ode {


/*-----------------------------------------------------------------------------
 *    
 *    @brief Abstract interface for Euler solver (implicit and explicit) 
 *    of (ODE) RHS problem 
 *    
 *    @Marco Ghiani Dec 2017, Glasgow UK
 *
 -----------------------------------------------------------------------------*/




template<typename Type = double>
class Euler :
                            public OdeSolver<Type> 
{
      
    public:  
      Euler(const rhsOdeProblem<Type> & that) noexcept :
                                                                        OdeSolver<Type>{that} 
                  {}
      
      virtual ~Euler() = default ;

      using OdeSolver<Type>::rhs;

      void solve(const std::string filename) override = 0  ;
      void solve() noexcept override  =0                   ;
//
//
   protected:
      Type uOld ;
      Type uNew ;
      Type  err ; 


};



  
  }//ode
 }//numeric
}//mg
# endif
