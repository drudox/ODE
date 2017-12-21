# ifndef __ADAMS_MOULTON_3STEP_4TH_ORDER_SOLVER_H__
# define __ADAMS_MOULTON_3STEP_4TH_ORDER_SOLVER_H__

# include "../AdamsMethods.H"
# include "../../../rhsOdeProblem.H"

namespace mg {
                namespace numeric {
                                     namespace ode {


/**-----------------------------------------------------------------------------------
 * @class AdamsMoulton4thSolver
 * @brief Perform implicit multi-step PREDICTOR - CORRECTOR method 
 * 4th order accuracy (4 step predictor - 3 step corrector ) solution of non stiff ODE  
 * 
 *
 * @author Marco Ghiani 
 * @date Dec 2017, Glasgow UK
 *
 *
 -------------------------------------------------------------------------------------*/




template<typename Type = double>
class AdamsMoulton4thSolver :
                                 public AdamsMethods<Type> 
{
      
    public:  
      AdamsMoulton4thSolver(const rhsOdeProblem<Type> & that) noexcept :
                                                                        AdamsMethods<Type>{that} 
                  {}
      
      virtual ~AdamsMoulton4thSolver() = default ;

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
      using AdamsMethods<Type>::k3 ;
      using AdamsMethods<Type>::k4 ;

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
inline void AdamsMoulton4thSolver<Type>::solve(const std::string filename) {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {  
         std::string mess = "Error opening file: " + filename + " in ADAMS-MOULTON 4th order solver" ;   
         throw std::runtime_error(mess.c_str());
      }
      else
      {
         std::cout << "Running Adams Bashforth (4step), CORRECTOR: Adams Moulton 4th order solver" << std::endl;
      
         t.at(0) = t0();
         u.at(0) = u0();   // initial Value 
         
         f << t.at(0) << ' ' << u.at(0) << std::endl; 
            
         // compute first point  (start-up the solver)  
         k1 = rhs.f(t.at(0)          , u.at(0));
         k2 = rhs.f(t.at(0)+ dt()/2. , u.at(0) + k1*dt()/2. );
         k3 = rhs.f(t.at(0)+ dt()/2. , u.at(0) + k2*dt()/2. );
         k4 = rhs.f(t.at(0)+ dt()    , u.at(0) + k3*dt() );
         
         u.at(1) =  u.at(0) + dt()/6.0 *(k1+2.*k2+2.*k3+k4);  // RungeKutta 4th order 
   
         //  first step 
         t.at(1) = t.at(0) + dt() ;
   
         f << t.at(1) << ' ' << u.at(1) << std::endl; // write to file first point 
 
         // compute second point  (start-up the solver)  
         k1 = rhs.f(t.at(1)          , u.at(1)              );
         k2 = rhs.f(t.at(1)+ dt()/2. , u.at(1) + k1*dt()/2. );
         k3 = rhs.f(t.at(1)+ dt()/2. , u.at(1) + k2*dt()/2. );
         k4 = rhs.f(t.at(1)+ dt()    , u.at(1) + k3*dt()    );
         
         u.at(2) =  u.at(1) + dt()/6.0 *(k1+2.*k2+2.*k3 + k4);  // RungeKutta 4th order order 
   
         //  second step 
         t.at(2) = t.at(1) + dt() ;
   
         f << t.at(2) << ' ' << u.at(2) << std::endl; // write to file first point 

         // compute third point  (start-up the solver)  
         k1 = rhs.f(t.at(2)          , u.at(2)              );
         k2 = rhs.f(t.at(2)+ dt()/2. , u.at(2) + k1*dt()/2. );
         k3 = rhs.f(t.at(2)+ dt()/2. , u.at(2) + k2*dt()/2. );
         k4 = rhs.f(t.at(2)+ dt()    , u.at(2) + k3*dt()    );
         
         u.at(3) =  u.at(2) + dt()/6.0 *(k1+2.*k2+2.*k3 + k4);  // RungeKutta 4th order order 
   
         //  third step 
         t.at(3) = t.at(2) + dt() ;
   
         f << t.at(3) << ' ' << u.at(3) << std::endl; // write to file third point 






         for(auto i=3; i < Ns ; i++ )
         {
            t.at(i+1) = t.at(i) + dt() ;

            // PREDICTOR
            //
            uPred = u.at(i) + dt()/24. *(55. * rhs.f(t.at(i), u.at(i)) 
                                        -59. * rhs.f(t.at(i-1), u.at(i-1)) 
                                        +37. * rhs.f(t.at(i-2), u.at(i-2)) 
                                        + 9. * rhs.f(t.at(i-3), u.at(i-3)) ) ;
            
            fPred = rhs.f(t.at(i+1),uPred);

            // CORRECTOR ADAMS MOULTON 
            
            error = 1.0 ;
            
            fCorrOld = fPred ;

          std::size_t iter = 0;  
              while(error >= pcToll)
              {
                 uCorrOld = u.at(i) + dt()/24.0 * ( 9. * fCorrOld 
                                                   +19. * rhs.f(t.at(i)  ,u.at(i)) 
                                                   -5. * rhs.f(t.at(i-1),u.at(i-1) )
                                                   +1. * rhs.f(t.at(i-2),u.at(i-2) ));
            
                 fCorrOld = rhs.f(t.at(i+1), uCorrOld) ;

                 uCorr    = u.at(i) + dt()/24.0 * (  9. * fCorrOld +
                                                   +19. * rhs.f(t.at(i)  , u.at(i)   )
                                                   -5.  * rhs.f(t.at(i-1), u.at(i-1) )  
                                                   +1.  * rhs.f(t.at(i-2), u.at(i-2) )  );

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
inline void AdamsMoulton4thSolver<Type>::solve() noexcept 
{

 std::cout << "Running Adams Bashforth (4step), CORRECTOR: Adams Moulton 4th order solver" << std::endl;
      
         t.at(0) = t0();
         u.at(0) = u0();   // initial Value 
         
         std::cout << t.at(0) << ' ' << u.at(0) << std::endl; 
            
         // compute first point  (start-up the solver)  
         k1 = rhs.f(t.at(0)          , u.at(0));
         k2 = rhs.f(t.at(0)+ dt()/2. , u.at(0) + k1*dt()/2. );
         k3 = rhs.f(t.at(0)+ dt()/2. , u.at(0) + k2*dt()/2. );
         k4 = rhs.f(t.at(0)+ dt()    , u.at(0) + k3*dt() );
         
         u.at(1) =  u.at(0) + dt()/6.0 *(k1+2.*k2+2.*k3+k4);  // RungeKutta 4th order 
   
         //  first step 
         t.at(1) = t.at(0) + dt() ;
   
         std::cout<< t.at(1) << ' ' << u.at(1) << std::endl; // write to file first point 
 
         // compute second point  (start-up the solver)  
         k1 = rhs.f(t.at(1)          , u.at(1)              );
         k2 = rhs.f(t.at(1)+ dt()/2. , u.at(1) + k1*dt()/2. );
         k3 = rhs.f(t.at(1)+ dt()/2. , u.at(1) + k2*dt()/2. );
         k4 = rhs.f(t.at(1)+ dt()    , u.at(1) + k3*dt()    );
         
         u.at(2) =  u.at(1) + dt()/6.0 *(k1+2.*k2+2.*k3 + k4);  // RungeKutta 4th order order 
   
         //  second step 
         t.at(2) = t.at(1) + dt() ;
   
         std::cout << t.at(2) << ' ' << u.at(2) << std::endl; // write to file first point 

         // compute third point  (start-up the solver)  
         k1 = rhs.f(t.at(2)          , u.at(2)              );
         k2 = rhs.f(t.at(2)+ dt()/2. , u.at(2) + k1*dt()/2. );
         k3 = rhs.f(t.at(2)+ dt()/2. , u.at(2) + k2*dt()/2. );
         k4 = rhs.f(t.at(2)+ dt()    , u.at(2) + k3*dt()    );
         
         u.at(3) =  u.at(2) + dt()/6.0 *(k1+2.*k2+2.*k3 + k4);  // RungeKutta 4th order order 
   
         //  third step 
         t.at(3) = t.at(2) + dt() ;
   
         std::cout << t.at(3) << ' ' << u.at(3) << std::endl; // write to file third point 






         for(auto i=3; i < Ns ; i++ )
         {
            t.at(i+1) = t.at(i) + dt() ;

            // PREDICTOR
            //
            uPred = u.at(i) + dt()/24. *(55. * rhs.f(t.at(i), u.at(i)) 
                                        -59. * rhs.f(t.at(i-1), u.at(i-1)) 
                                        +37. * rhs.f(t.at(i-2), u.at(i-2)) 
                                        + 9. * rhs.f(t.at(i-3), u.at(i-3)) ) ;
            
            fPred = rhs.f(t.at(i+1),uPred);

            // CORRECTOR ADAMS MOULTON 
            
            error = 1.0 ;
            
            fCorrOld = fPred ;

          std::size_t iter = 0;  
              while(error >= pcToll)
              {
                 uCorrOld = u.at(i) + dt()/24.0 * ( 9. * fCorrOld 
                                                   +19. * rhs.f(t.at(i)  ,u.at(i)) 
                                                   -5. * rhs.f(t.at(i-1),u.at(i-1) )
                                                   +1. * rhs.f(t.at(i-2),u.at(i-2) ));
            
                 fCorrOld = rhs.f(t.at(i+1), uCorrOld) ;

                 uCorr    = u.at(i) + dt()/24.0 * (  9. * fCorrOld +
                                                   +19. * rhs.f(t.at(i)  , u.at(i)   )
                                                   -5.  * rhs.f(t.at(i-1), u.at(i-1) )  
                                                   +1.  * rhs.f(t.at(i-2), u.at(i-2) )  );

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
