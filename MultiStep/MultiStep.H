# ifndef __MULTISTEP_H__
# define __MULTISTEP_H__ 

# include "../rhsOdeProblem.H"
# include "../OdeSolver.H"


namespace mg { 
                namespace numeric {
                                    namespace ode {


/*-------------------------------------------------------------------------------
 *    
 *    Base class for all the multistep solvers of IVP :
 *    dy/dt = f(y,t)
 *
 *    @Marco Ghiani Dec 2017, Glasgow UK
 *
 ------------------------------------------------------------------------------*/


template <typename Type = double> 
class MultiStep :      
                        public OdeSolver<Type>
{
    
   public: 
      
      MultiStep(const rhsOdeProblem<Type>& that ) noexcept : 
                                                               OdeSolver<Type>{that} 
                          {}                                
      
      virtual ~MultiStep() = default ;
      
      using OdeSolver<Type>::rhs;
      

      virtual void solve(const std::string filename) override = 0;
      virtual void solve() noexcept  override                 = 0;


   protected:
     
      std::vector<Type> up1 ;
      std::vector<Type> um1 ;
      std::vector<Type> um2 ;

};



  }//ode 
 }//numeric
}//mg
# endif 
