# ifndef __LEAPFROG_SOLVER_H__
# define __LEAPFROG_SOLVER_H__

# include "MultiStep.H"
# include "../rhsOdeProblem.H"

namespace mg {
                namespace numeric {
                                     namespace ode {


/*-----------------------------------------------------------------------------
 *    Compute semi-Implicit solution  of a given (ODE) RHS problem using 
 *    
 *    - LEAP-FROG (mid-point) scheme 
 *    
 *    @Marco Ghiani Nov 2017, Glasgow UK
 *
 -----------------------------------------------------------------------------*/




template<typename Type = double>
class LeapFrogSolver :                                      // LEAP-FROG SOLVER
                            public MultiStep<Type> 
{
      
    public:  
      LeapFrogSolver(const rhsOdeProblem<Type> & that) noexcept :
                                                                    MultiStep<Type>{that} 
                  {}
      
      virtual ~LeapFrogSolver() = default ;

      using OdeSolver<Type>::rhs;

      void solve(const std::string filename) override final ;
      void solve() noexcept override final                  ;
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
inline void LeapFrogSolver<Type>::solve(const std::string filename) {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {     
         throw std::runtime_error(">> Error opening File in FWD euler solver <<" );
      }
      else
      {
         std::cout << "Running LeapFrog (Leap-Frog) Solver" << std::endl;
      
         t.at(0) = t0();
         u.at(0) = u0();
         
         f << t.at(0) << ' ' << u.at(0) << std::endl; 
        
            
         //    
         Type k1 = rhs.f(t.at(0) , u.at(0) );
         Type k2 = rhs.f(t.at(0)+ dt()/2 , u.at(0) + k1*dt()/2 );
   
         // initiation first point 
         u.at(1) =  u.at(0) + dt() * k2;  // Rk 2nd order PREDICTOR
         //u.at(1) = u.at(0) + dt() * rhs.f(t.at(0) , u.at(0)) ; // FWD EULER PREDICTOR    
         t.at(1) = t.at(0) + dt() ;
         
         f << t.at(1) << ' ' << u.at(1) << std::endl;     

         for(auto i=1; i < Ns ; i++ )
         {
            t.at(i+1) = t.at(i) + dt() ;
            u.at(i+1) = u.at(i-1) + 2*dt() * rhs.f(t.at(i), u.at(i)) ;  // leap-frog 
            f << t.at(i+1) << ' ' << u.at(i+1) << std::endl ;
         } 
         std::cout << "... Done " << std::endl;  
      
         f.close();
      } 
}


template<typename Type>
inline void LeapFrogSolver<Type>::solve() noexcept 
{
     std::cout << "Running LeapFrog (Leap-Frog) Solver" << std::endl;
       t.at(0) = t0();
       u.at(0) = u0();
         
       std::cout << t.at(0) << ' ' << u.at(0) << std::endl; 
       
       // initiation first point 
       u.at(1) = u.at(0) + dt() * rhs.f(t.at(0) , u.at(0)) ; // FWD EULER PREDICTOR    
       t.at(1) = t.at(0) + dt() ;
  
       for(auto i=1; i < Ns ; i++ )
       {
          t.at(i+1) = t.at(i) + dt() ;
          u.at(i+1) = u.at(i-1) + 2*dt() * rhs.f(t.at(i), u.at(i)) ;  // leap frog
          std::cout << t.at(i) << ' ' << u.at(i) << std::endl ;
       } 
       std::cout << "... Done " << std::endl;  
}
  
  }//ode
 }//numeric
}//mg
# endif
