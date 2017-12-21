# ifndef __ADAMS_MOULTON_2STEP_3TH_ORDER_SOLVER_H__
# define __ADAMS_MOULTON_2STEP_3TH_ORDER_SOLVER_H__

# include "../AdamsMethods.H"
# include "../../../rhsOdeProblem.H"

namespace mg {
                namespace numeric {
                                     namespace ode {


/**-----------------------------------------------------------------------------------
 * @class AdamsMoulton3thSolver
 * @brief 
 * Perform implicit multi-step 3th order accuracy (2 step) solution of a given (ODE).  
 * 
 *
 * @author Marco Ghiani 
 * @date Dec 2017, Glasgow UK
 *
 *
 -------------------------------------------------------------------------------------*/




template<typename Type = double>
class AdamsMoulton3thSolver :
                                 public AdamsMethods<Type> 
{
      
    public:  
      AdamsMoulton3thSolver(const rhsOdeProblem<Type> & that) noexcept :
                                                                        AdamsMethods<Type>{that} 
                  {}
      
      virtual ~AdamsMoulton3thSolver() = default ;

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

      using AdamsMethods<Type>::k1 ;
      using AdamsMethods<Type>::k2 ;

      using AdamsMethods<Type>::uPred ;
      using AdamsMethods<Type>::uCorr ;

      using AdamsMethods<Type>::uCorrOld ;
      using AdamsMethods<Type>::fPred ;
      using AdamsMethods<Type>::fCorr ;
      using AdamsMethods<Type>::fCorrOld ;
      using AdamsMethods<Type>::error ;
      
      using AdamsMethods<Type>::pcToll ;


};

//------------------  Implementation (to be put into .cpp file)   -----------------  //


template<typename Type>
inline void AdamsMoulton3thSolver<Type>::solve(const std::string filename) {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {  
         std::string mess = "Error opening file: " + filename + " in Solver ADAMS-MOULTON 3th order" ;   
         throw std::runtime_error(mess.c_str());
      }
      else
      {
         std::cout << "Running Adams Bashforth (3step), CORRECTOR: Adams Moulton 3th order solver" << std::endl;
      
         t.at(0) = t0();
         u.at(0) = u0();   // initial Value 
         
         f << t.at(0) << ' ' << u.at(0) << std::endl; 
            
         // compute first point  (start-up the solver)  
         k1 = rhs.f(t.at(0) , u.at(0));
         k2 = rhs.f(t.at(0)+ dt() , u.at(0) + k1*dt() );
         
         u.at(1) =  u.at(0) + dt()/2 *(k1+k2);  // RungeKutta 2nd order 
   
         //  first point 
         t.at(1) = t.at(0) + dt() ;
   
         f << t.at(1) << ' ' << u.at(1) << std::endl; // write to file first point 
 
         // compute second point  (start-up the solver)  
         k1 = rhs.f(t.at(1) , u.at(1));
         k2 = rhs.f(t.at(1)+ dt() , u.at(1) + k1*dt() );
         
         u.at(2) =  u.at(1) + dt()/2 *(k1+k2);  // RungeKutta 2nd order 
   
         //  second point 
         t.at(2) = t.at(1) + dt() ;
   
         f << t.at(2) << ' ' << u.at(2) << std::endl; // write to file first point 


         for(auto i=2; i < Ns ; i++ )
         {
            t.at(i+1) = t.at(i) + dt() ;

            // PREDICTOR
            //
            uPred = u.at(i) + dt()/12. *(23. * rhs.f(t.at(i), u.at(i)) 
                                        -16. * rhs.f(t.at(i-1), u.at(i-1)) 
                                        + 5. * rhs.f(t.at(i-2), u.at(i-2)) ) ;
            
            fPred = rhs.f(t.at(i+1),uPred);

            // CORRECTOR ADAMS MOULTON 
            
            error = 1.0 ;
            
            fCorrOld = fPred ;

          std::size_t iter = 0;  
              while(error >= pcToll)
              {
                 uCorrOld = u.at(i) + dt()/12.0 * ( 5. * fCorrOld 
                                                   +8. * rhs.f(t.at(i)  ,u.at(i)) 
                                                   -1. * rhs.f(t.at(i-1),u.at(i-1) ));
            
                 fCorrOld = rhs.f(t.at(i+1), uCorrOld) ;

                 uCorr    = u.at(i) + dt()/12.0 * ( 5. * fCorrOld 
                                                   +8. * rhs.f(t.at(i)  , u.at(i)   )
                                                   -1. * rhs.f(t.at(i-1), u.at(i-1) )  );

                 fCorr    = rhs.f(t.at(i+1) , uCorr );
            
                 error    = fabs(uCorr-uCorrOld);

                 fCorrOld = fCorr ;
            
                 iter++;
              } 
            
            u.at(i+1) = uCorr ;

            f << t.at(i+1) << ' ' << u.at(i+1) << std::endl ;
         
         }
         std::cout << "... Done " << std::endl;  
      
         f.close();
      } 
}


template<typename Type>
inline void AdamsMoulton3thSolver<Type>::solve() noexcept 
{
  std::cout << "Running Adams Bashforth (3step), CORRECTOR: Adams Moulton 3 step Solver" << std::endl;
      
         t.at(0) = t0();
         u.at(0) = u0();   // initial Value 
         
         std::cout << t.at(0) << ' ' << u.at(0) << std::endl; 
            
         // compute first point  (start-up the solver)  
         k1 = rhs.f(t.at(0) , u.at(0));
         k2 = rhs.f(t.at(0)+ dt() , u.at(0) + k1*dt() );
         
         u.at(1) =  u.at(0) + dt()/2 *(k1+k2);  // RungeKutta 2nd order 
   
         //  first point 
         t.at(1) = t.at(0) + dt() ;
   
         std::cout << t.at(1) << ' ' << u.at(1) << std::endl; // write to file first point 
 
         // compute second point  (start-up the solver)  
         k1 = rhs.f(t.at(1) , u.at(1));
         k2 = rhs.f(t.at(1)+ dt() , u.at(1) + k1*dt() );
         
         u.at(2) =  u.at(1) + dt()/2 *(k1+k2);  // RungeKutta 2nd order 
   
         //  second point 
         t.at(2) = t.at(1) + dt() ;
   
         std::cout << t.at(2) << ' ' << u.at(2) << std::endl; // write to file first point 


         for(auto i=2; i < Ns ; i++ )
         {
            t.at(i+1) = t.at(i) + dt() ;

            // PREDICTOR
            //
            uPred = u.at(i) + dt()/12. *(23. * rhs.f(t.at(i), u.at(i)) 
                                        -16. * rhs.f(t.at(i-1), u.at(i-1)) 
                                        + 5. * rhs.f(t.at(i-2), u.at(i-2)) ) ;
            
            fPred = rhs.f(t.at(i+1),uPred);

            // CORRECTOR ADAMS MOULTON 
            
            error = 1.0 ;
            
            fCorrOld = fPred ;

          std::size_t iter = 0;  
              while(error >= pcToll)
              {
                 uCorrOld = u.at(i) + dt()/12.0 * ( 5. * fCorrOld 
                                                   +8. * rhs.f(t.at(i)  ,u.at(i)) 
                                                   -1. * rhs.f(t.at(i-1),u.at(i-1) ));
            
                 fCorrOld = rhs.f(t.at(i+1), uCorrOld) ;

                 uCorr    = u.at(i) + dt()/12.0 * ( 5. * fCorrOld +
                                                    8. * rhs.f(t.at(i)  , u.at(i)   )
                                                   -1. * rhs.f(t.at(i-1), u.at(i-1) )  );

                 fCorr    = rhs.f(t.at(i+1) , uCorr );
            
                 error    = fabs(uCorr-uCorrOld);

                 fCorrOld = fCorr ;
            
                 iter++;
              } 
            
            u.at(i+1) = uCorr ;

            std::cout << t.at(i+1) << ' ' << u.at(i+1) << std::endl ;
         
         }
     
        std::cout << "... Done " << std::endl;  
}
  
  }//ode
 }//numeric
}//mg
# endif
