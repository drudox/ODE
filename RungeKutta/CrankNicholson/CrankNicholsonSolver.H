# ifndef __CRANK_NICHOLSON_SOLVER_H__
# define __CRANK_NICHOLSON_SOLVER_H__

# include "../RungeKutta.H"
# include "../../rhsOdeProblem.H"

namespace mg {
                namespace numeric {
                                     namespace ode {


/*-----------------------------------------------------------------------------
 *    
 *     Compute Crank-Nicholson Implicit scheme for IVP 
 *    
 *     using >> Predictor-Corrector method
 *
 *    
 *    @author Marco Ghiani Dec 2017, Glasgow UK
 *
 -----------------------------------------------------------------------------*/




template<typename Type = double>
class CrankNicholsonSolver :
                            public RungeKutta<Type> 
{
      
    public:  
      CrankNicholsonSolver(const rhsOdeProblem<Type> & that) noexcept :
                                                                        RungeKutta<Type>{that} 
                  {}
      
      virtual ~CrankNicholsonSolver() = default ;

      using OdeSolver<Type>::rhs;

      void solve(const std::string filename) override final;
      void solve() noexcept override final                 ;
//
//
  private:

      using OdeSolver<Type>::t  ;
      using OdeSolver<Type>::u  ;
      
      using RungeKutta<Type>::uc ;

      using OdeSolver<Type>::dt ; 
      using OdeSolver<Type>::t0 ;
      using OdeSolver<Type>::tf ;
      using OdeSolver<Type>::u0 ;
      
      using OdeSolver<Type>::Ns ;

};

//------------------  Implementation (to be put into .cpp file)   -----------------  //


template<typename Type>
inline void CrankNicholsonSolver<Type>::solve(const std::string filename) {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {     
         throw std::runtime_error(">> Error opening File in FWD euler solver <<" );
      }
      else
      {
         std::cout << "Running CrankNicholson Solver" << std::endl;
      
         t.at(0)  = t0();
         u.at(0)  = u0();
         uc.at(0) = u0();
            

         f << t.at(0) << ' ' << uc.at(0) << std::endl; 
         
         for(auto i=1; i <= Ns ; i++ )
         {
            t.at(i) = t.at(i-1) + dt() ;
            // Predictor step (Exp Euler) 
            u.at(i) = u.at(i-1) + dt() *rhs.f(t.at(i-1), u.at(i-1)) ;  
            // Corrector step : Crank Nicolson 
            uc.at(i) = uc.at(i-1) + dt()/2 * ( rhs.f(t.at(i-1),uc.at(i-1)) + rhs.f(t.at(i), u.at(i)) ) ;

            f << t.at(i) << ' ' << uc.at(i) << std::endl ;
         
         } 
         std::cout << "... Done " << std::endl;  
      
         f.close();
      } 
}


template<typename Type>
inline void CrankNicholsonSolver<Type>::solve() noexcept 
{
       std::cout << "Running CrankNicholson Solver" << std::endl;
      
         t.at(0)  = t0();
         u.at(0)  = u0();
         uc.at(0) = u0();
         std::cout << t.at(0) << ' ' << uc.at(0) << std::endl; 
         
         for(auto i=1; i <= Ns ; i++ )
         {
            t.at(i) = t.at(i-1) + dt() ;
            // Predictor step (Exp Euler) 
            u.at(i) = u.at(i-1) + dt() *rhs.f(t.at(i-1), u.at(i-1)) ;  
            // Corrector step : Crank Nicolson 
            uc.at(i) = uc.at(i-1) + dt()/2 * ( rhs.f(t.at(i-1),uc.at(i-1)) + rhs.f(t.at(i), u.at(i)) ) ;

            std::cout << t.at(i) << ' ' << uc.at(i) << std::endl ;
         
         } 
         std::cout << "... Done " << std::endl;  
      
 
}
  
  }//ode
 }//numeric
}//mg
# endif
