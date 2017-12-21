# ifndef __ODE_SOLVER_H__
# define __ODE_SOLVER_H__

# include "abstractOdeSolver.H"
# include <cmath>
# include <iostream>
# include <cstdlib>
# include "rhsOdeProblem.H"
# include <vector>
# include <string>
# include "rhsOdeProblem.H"

namespace mg { 
                namespace numeric {
                                     namespace ode {

template< typename T>
class rhsOdeProblem;


/*------------------------------------------------------------*
 *    
 *    Abstract OdeSolver Class - double or single precision  
 *    
 *    @Marco Ghiani , Oct 2017 Glasgow
 *
 -------------------------------------------------------------*/


template <typename Type>
class OdeSolver  : 
                    public AbstractODESolver<Type> 
{
  
//--
//
  public:

    OdeSolver(const rhsOdeProblem<Type>& that) noexcept : rhs{that} 
    {
     
     setStepSize    () ; 
     setInitialTime () ;
     setFinalTime   () ;
     setInitialValue() ;
     setSize() ;
    }

    virtual ~OdeSolver() = 0;

    rhsOdeProblem<Type> rhs ;
    

     virtual void setStepSize    () override { stepSize     = rhs.dt(); }
     virtual void setInitialTime () override { initialTime  = rhs.t0(); }
     virtual void setFinalTime   () override { finalTime    = rhs.tf(); }
     virtual void setInitialValue() override { initialValue = rhs.u0(); }

     virtual void solve(const std::string filename)     = 0;
     virtual void solve() noexcept                      = 0;

     
     auto setRhs(const rhsOdeProblem<Type>& that ) { this->rhs = that ; }  
     
     virtual Type getStepSize()    const override { return stepSize     ;}
     virtual Type getInitialTime() const override { return initialTime  ;}
     virtual Type getFinalTime()   const override { return finalTime    ;}
     virtual Type getInitialValue()const override { return initialValue ;}
      
     virtual Type dt() const noexcept { return stepSize     ;}
     virtual Type t0() const noexcept { return initialTime  ;}
     virtual Type tf() const noexcept { return finalTime    ;}
     virtual Type u0() const noexcept { return initialValue ;}
     
     virtual void setSize() noexcept ;

     protected:
      
      Type stepSize;
      Type initialTime;
      Type finalTime;
      Type initialValue;
      
      const int Ns = (rhs.tf() - rhs.t0())/rhs.dt() ;
      
      constexpr static Type toll = 1e-12 ;
      
      std::vector<Type> t ;
      std::vector<Type> u ;
      std::vector<Type> uc ;
      
};

template<typename Type>
OdeSolver<Type>::~OdeSolver() = default ;

template<typename Type>
void OdeSolver<Type>::setSize() noexcept
{
  t.resize(Ns+1) ;  
  u.resize(Ns+1) ;    
  uc.resize(Ns+1) ;    
}


  }//ode
 }//numeric
}//mg
# endif
