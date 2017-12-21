# ifndef __ADAMS_METHODS_H__
# define __ADAMS_METHODS_H__ 

# include "../../rhsOdeProblem.H"
# include "../MultiStep.H"


namespace mg { 
                namespace numeric {
                                    namespace ode {


/*-------------------------------------------------------------------------------
 *    
 *    @brief Base class for all the Adams (Bashforth - Moulton) multistep solver
 *    for ODE problems :
 *    dy/dt = f(y,t)
 *
 *    @author Marco Ghiani Dec 2017, Glasgow UK
 *
 ------------------------------------------------------------------------------*/


template <typename Type = double> 
class AdamsMethods :      
                        public MultiStep<Type>
{
    
   public: 
      
      AdamsMethods(const rhsOdeProblem<Type>& that ) noexcept : 
                                                               MultiStep<Type>{that} 
                          {}                                
      
      virtual ~AdamsMethods() = default ;
      
      using OdeSolver<Type>::rhs;
      

      virtual void solve(const std::string filename) override = 0;
      virtual void solve() noexcept  override                 = 0;

      //void reSize() noexcept ;

   protected:
     
     /*
      std::vector<Type> up1 ;
      std::vector<Type> um1 ;
      std::vector<Type> um2 ;
      */

     Type k1 ; 
     Type k2 ; 
     Type k3 ; 
     Type k4 ; 
     Type k5 ; 
      
     Type uPred;
     Type uCorr;
     Type uCorrOld ;
     Type fPred ;
     Type fCorr ;
     Type fCorrOld ;
     Type error ;
     
     constexpr static Type pcToll = 1e-10;

};



  }//ode 
 }//numeric
}//mg
# endif 
