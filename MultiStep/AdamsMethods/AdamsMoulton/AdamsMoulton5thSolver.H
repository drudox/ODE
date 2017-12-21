# ifndef __ADAMS_MOULTON_4STEP_5TH_ORDER_SOLVER_H__
# define __ADAMS_MOULTON_4STEP_5TH_ORDER_SOLVER_H__

# include "../AdamsMethods.H"
# include "../../../rhsOdeProblem.H"

namespace mg {
                namespace numeric {
                                     namespace ode {


/**--------------------------------------------------------------------------------------------
 * @class AdamsMoulton5thSolver
 * @brief Perform implicit multi-step 5th order accuracy using PREDICTOR - CORRECTOR method
 * (5 step predictor - 4 step implicit corrector)
 * Compute the solution of a given non-stiff ODE.  
 * 
 *
 * @author Marco Ghiani 
 * @date Dec 2017, Glasgow UK
 *
 *
 -------------------------------------------------------------------------------------*/




template<typename Type = double>
class AdamsMoulton5thSolver :
                                 public AdamsMethods<Type> 
{
      
    public:  
      AdamsMoulton5thSolver(const rhsOdeProblem<Type> & that) noexcept :
                                                                        AdamsMethods<Type>{that} 
                  {}
      
      virtual ~AdamsMoulton5thSolver() = default ;

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
      using AdamsMethods<Type>::k5 ;

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
inline void AdamsMoulton5thSolver<Type>::solve(const std::string filename) {
      
      std::ofstream f(filename, std::ios::out );
      
      if(!f)
      {  
         std::string mess = "Error opening file: " + filename + " in ADAMS-MOULTON 5th order solver" ;   
         throw std::runtime_error(mess.c_str());
      }
      else
      {
         std::cout << "Running Adams Bashforth (5step), CORRECTOR: Adams Moulton 5th order solver" << std::endl;
      
         t.at(0) = t0();
         u.at(0) = u0();   // initial Value 
         
         f << t.at(0) << ' ' << u.at(0) << std::endl; 
            
         // compute first point  (start-up the solver)  
         k1 = dt()*rhs.f(t.at(0)          , u.at(0)                         );
         k2 = dt()*rhs.f(t.at(0)+ dt()/3. , u.at(0) + k1/3.                 );
         k3 = dt()*rhs.f(t.at(0)+ dt()/3. , u.at(0) + 1./6.*(k1+k2)         );
         k4 = dt()*rhs.f(t.at(0)+ dt()/2. , u.at(0) + 1./8.*(k1+3.*k3)      );
         k5 = dt()*rhs.f(t.at(0)+ dt()    , u.at(0) + 1./2.*(k1-3.*k3+4.*k4));
         
         u.at(1) =  u.at(0) + 1./6.0 *(k1+4.*k4+1.*k5);  // Runge-Kutta-Merson 5th order 
   
         //  first step 
         t.at(1) = t.at(0) + dt() ;
   
         f << t.at(1) << ' ' << u.at(1) << std::endl; // write to file first point 
 
         // compute second point  (start-up the solver)  
         k1 = dt()*rhs.f(t.at(1)          , u.at(1)                         );
         k2 = dt()*rhs.f(t.at(1)+ dt()/3. , u.at(1) + k1/3.                 );
         k3 = dt()*rhs.f(t.at(1)+ dt()/3. , u.at(1) + 1./6.*(k1+k2)         );
         k4 = dt()*rhs.f(t.at(1)+ dt()/2. , u.at(1) + 1./8.*(k1+3.*k3)      );
         k5 = dt()*rhs.f(t.at(1)+ dt()    , u.at(1) + 1./2.*(k1-3.*k3+4.*k4));
         
         u.at(2) =  u.at(1) + 1./6.0 *(k1+4.*k4+1.*k5);  // Runge-Kutta-Merson 5th order 
   
         //  second step 
         t.at(2) = t.at(1) + dt() ;
   
         f << t.at(2) << ' ' << u.at(2) << std::endl; // write to file first point 

         // compute third point  (start-up the solver)  
         k1 = dt()*rhs.f(t.at(2)          , u.at(2)                         );
         k2 = dt()*rhs.f(t.at(2)+ dt()/3. , u.at(2) + k1/3.                 );
         k3 = dt()*rhs.f(t.at(2)+ dt()/3. , u.at(2) + 1./6.*(k1+k2)         );
         k4 = dt()*rhs.f(t.at(2)+ dt()/2. , u.at(2) + 1./8.*(k1+3.*k3)      );
         k5 = dt()*rhs.f(t.at(2)+ dt()    , u.at(2) + 1./2.*(k1-3.*k3+4.*k4));
         
         u.at(3) =  u.at(2) + 1./6.0 *(k1+4.*k4+1.*k5);  // Runge-Kutta-Merson 5th order 
   
        //  third step 
         t.at(3) = t.at(2) + dt() ;
   
         f << t.at(3) << ' ' << u.at(3) << std::endl; // write to file third point 

         // compute fourth point  (start-up the solver)  
         k1 = dt()*rhs.f(t.at(3)          , u.at(3)                         );
         k2 = dt()*rhs.f(t.at(3)+ dt()/3. , u.at(3) + k1/3.                 );
         k3 = dt()*rhs.f(t.at(3)+ dt()/3. , u.at(3) + 1./6.*(k1+k2)         );
         k4 = dt()*rhs.f(t.at(3)+ dt()/2. , u.at(3) + 1./8.*(k1+3.*k3)      );
         k5 = dt()*rhs.f(t.at(3)+ dt()    , u.at(3) + 1./2.*(k1-3.*k3+4.*k4));
         
         u.at(4) =  u.at(3) + 1./6.0 *(k1+4.*k4+1.*k5);  // Runge-Kutta-Merson 5th order 
   
        //  fourth step 
         t.at(4) = t.at(3) + dt() ;
   
         f << t.at(4) << ' ' << u.at(4) << std::endl; // write to file third point 






         for(auto i=4; i < Ns ; i++ )
         {
            t.at(i+1) = t.at(i) + dt() ;

            // PREDICTOR
            //
            uPred = u.at(i) + dt()     *(1901.0/720.0 * rhs.f(t.at(i), u.at(i)) 
                                        -1387.0/360.0 * rhs.f(t.at(i-1), u.at(i-1)) 
                                        +109.0/30.0   * rhs.f(t.at(i-2), u.at(i-2)) 
                                        -637.0/360.0  * rhs.f(t.at(i-3), u.at(i-3)) 
                                        +251.0/720.0  * rhs.f(t.at(i-4), u.at(i-4)) ) ;
            
            fPred = rhs.f(t.at(i+1),uPred);

            // CORRECTOR ADAMS MOULTON 
            
            error = 1.0 ;
            
            fCorrOld = fPred ;

          std::size_t iter = 0;  
              while(error >= pcToll)
              {
                 uCorrOld = u.at(i) + dt()/720.0 * ( 251. * fCorrOld 
                                                    +646. * rhs.f(t.at(i)  ,u.at(i)) 
                                                    -264. * rhs.f(t.at(i-1),u.at(i-1) )
                                                    +106. * rhs.f(t.at(i-2),u.at(i-2) )
                                                    - 19. * rhs.f(t.at(i-3),u.at(i-3) ));
            
                 fCorrOld = rhs.f(t.at(i+1), uCorrOld) ;

                 uCorr    = u.at(i) + dt()/720.0 * ( 251. * fCorrOld 
                                                    +646. * rhs.f(t.at(i)  ,u.at(i)) 
                                                    -264. * rhs.f(t.at(i-1),u.at(i-1) )
                                                    +106. * rhs.f(t.at(i-2),u.at(i-2) )
                                                    - 19. * rhs.f(t.at(i-3),u.at(i-3) ));
 
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
inline void AdamsMoulton5thSolver<Type>::solve() noexcept 
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
