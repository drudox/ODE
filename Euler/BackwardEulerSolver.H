# ifndef __BACKWARD_EULER_SOLVER_H__
# define __BACKWARD_EULER_SOLVER_H__

# include "Euler.H"
# include "../rhsOdeProblem.H"

namespace mg { 
               namespace numeric {
                                    namespace ode {


/*----------------------------------------------------------------------*
 *    
 *    Compute Implicit Euler. Solve a given (ODE) RHS problem :
 *    y' = f(t,y); 
 *
 *    @author Marco Ghiani , 
 *    @date   Dec 2017 
 *    @place  Glasgow UK
 *
 -----------------------------------------------------------------------*/



template <typename Type = double>
class BackwardEulerSolver : 
                              public Euler<Type>
{
      
   public:  
      BackwardEulerSolver(const rhsOdeProblem<Type> & that) :
                                                               Euler<Type>{that} 
                  {}

      virtual ~BackwardEulerSolver() = default ;
      
      using OdeSolver<Type>::rhs;

      void solve(std::string filename) override final ;
      void solve()  noexcept override final ;
   
   private:
      
      using OdeSolver<Type>::t  ;
      using OdeSolver<Type>::u  ;
      using OdeSolver<Type>::dt ; 
      using OdeSolver<Type>::t0 ;
      using OdeSolver<Type>::tf ;
      using OdeSolver<Type>::u0 ;
      
      using OdeSolver<Type>::Ns ;
      
      using OdeSolver<Type>::toll ;
      
      using Euler<Type>::uOld;
      using Euler<Type>::uNew;
      using Euler<Type>::err ;
};

//------------------  Implementation (to be put into .cpp file) -------------------- //


template <typename Type>
inline void BackwardEulerSolver<Type>::solve(std::string filename)  {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {     
         std::string mess = "Error opening file " + filename + "in BackwardEulerSolver" ;  
         throw std::runtime_error(mess.c_str());         
      }
      else
      {
         std::cout << "Running BackwardEuler Solver" << std::endl;
      
         t.at(0) = t0();
         u.at(0) = u0();
         
         f << t.at(0) << ' ' << u.at(0) << std::endl; 
         
            
         auto max_it = 100;
         for(auto i=1; i <= Ns ; i++ )
         {
            t.at(i) = t.at(i-1) + dt() ;
            
            uOld = u.at(i-1) + dt() *rhs.f(t.at(i), u.at(i-1)) ;  
            
              err = 1.0;
            
            auto j = 1;

            while(err > toll ) //&& j < max_it)
            {
               uNew = uOld - ( uOld - (u.at(i-1) + dt() * rhs.f(t.at(i),uOld) ) ) /   
                             ( 1- dt() * rhs.dfdt(t.at(i),uOld) ) ;
               
               err = fabs(uNew - uOld );
               //std::cout << err << ' ' << j << std::endl;
               uOld = uNew;
               j++ ;
            }      
            
            u.at(i) = uNew ;
            
            f << t.at(i) << ' ' << u.at(i) << std::endl ;
         } 
         std::cout << "... Done " << std::endl;  
      f.close();
      }
}


template <typename Type>
inline void BackwardEulerSolver<Type>::solve() noexcept 
{
     std::cout << "Running BackwardEuler Solver" << std::endl;
 
         t.at(0) = t0();
         u.at(0) = u0();
         
         std::cout << t.at(0) << ' ' << u.at(0) << std::endl; 
         

         for(auto i=1; i <= Ns ; i++ )
         {
            t.at(i) = t.at(i-1) + dt() ;
            
            uOld = u.at(i-1) + dt() *rhs.f(t.at(i), u.at(i-1)) ;  
            
              err = 1.0;
            
              auto j = 1;

            while(err > toll)
            {
               uNew = uOld - ( uOld - (u.at(i-1) + dt() * rhs.f(t.at(i),uOld) ) ) /   
                             ( 1- dt() * rhs.dfdt(t.at(i),uOld) ) ;
               err = fabs(uNew - uOld );
               uOld = uNew;
               j++ ;
            }      
            
            u.at(i) = uNew ;
            
           std::cout << t.at(i) << ' ' << u.at(i) << std::endl ;
         } 
         std::cout << "... Done " << std::endl;  
}

  }//ode
 }//numeric
}//mg 
# endif
