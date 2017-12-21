# ifndef __RHS_ODE_PROBLEM_H__
# define __RHS_ODE_PROBLEM_H__

# include <string>
# include <functional>
# include <fstream>

namespace mg {
               namespace numeric {
                                    namespace ode {

/*-----------------------------------------------------------------------
 *   @brief Class RightHandSide - ODE problem,
 *    
 *    --> Type : single precision (float) , double precision (double)
 *
 *    dy/dt = RHS 
 *
 *    @ Marco Ghiani  Oct 2017 Glasgow UK
 ------------------------------------------------------------------------*/


template <typename Type = double>
class rhsOdeProblem {
  
//
//
//--  
   public:
      
      rhsOdeProblem (const std::function<const Type(const Type,const Type)> numfun ,
                     const std::function<const Type(const Type,const Type)> exactfun ,
                     const Type,const Type,const Type,const Type,const std::string ) noexcept ;
       
  
      rhsOdeProblem (const std::function<const Type(const Type,const Type)> numfun ,
                     const Type, const Type, const Type, const Type  ) noexcept ;
 
      
      virtual ~rhsOdeProblem() = default ;
      
      rhsOdeProblem(const rhsOdeProblem &) = default ;
      rhsOdeProblem(rhsOdeProblem&& ) = default ;
      rhsOdeProblem& operator=(const rhsOdeProblem&) = default ;
      rhsOdeProblem& operator=(rhsOdeProblem&& ) = default ;


      std::function<const Type(const Type,const Type)> numericalFunction ;
      
      std::function<const Type(const Type,const Type)> analiticalFunction ;

      Type    f(Type t, Type u) const noexcept { return numericalFunction(t,u); }
      Type dfdt(Type t, Type u) const noexcept { return (f(t,u+eps)-f(t,u))/eps ; }

      auto setRhs  (std::function<Type(Type,Type)> numfun) noexcept { numericalFunction = numfun; } 
      
      auto setExact(std::function<Type(Type,Type)> exactfun) noexcept {analiticalFunction = exactfun;}
    
      auto solveExact() noexcept ;

      const Type t0() const noexcept { return _t0 ;} 
      const Type tf() const noexcept { return _tf ;}
      const Type dt() const noexcept { return _dt ;}
      const Type u0() const noexcept { return _u0 ;}
      const std::string fname () const noexcept { return filename ;}
      

//---
   private:
     
     Type _t0 ;
     Type _tf ;
     Type _dt ;
     Type _u0 ;

     std::string filename ; 

     constexpr static Type eps  = 1e-12  ; 
      

};

/*
 *    Implementation 
 */ 

template <typename Type>
rhsOdeProblem<Type>::rhsOdeProblem ( const std::function<const Type(const Type,const Type)> numfun ,
                                     const std::function<const Type(const Type,const Type)> exactfun ,
                                     const Type Ti,const Type Tf,const Type Dt,const Type U0,
                                     const std::string fname 
                                   ) 
                                    noexcept : numericalFunction{numfun} ,
                                            analiticalFunction{exactfun} ,
                                                                 _t0{Ti} ,
                                                                 _tf{Tf} ,
                                                                 _dt{Dt} ,
                                                                 _u0{U0} ,
                                                         filename{fname} 
{
      solveExact();
}

template<typename Type>
rhsOdeProblem<Type>::rhsOdeProblem ( const std::function<const Type(const Type,const Type)> numfun ,
                                     const Type Ti,const Type Tf,const Type Dt,const Type U0              
                                   ) 
                                        noexcept : numericalFunction{numfun} , 
                                                                     _t0{Ti} ,
                                                                     _tf{Tf} ,
                                                                     _dt{Dt} ,
                                                                     _u0{U0}  
                   {}                                   


//- if exist ( and gives ) compute the 
//     numerical-exact solution 
//
template<typename Type>
auto rhsOdeProblem<Type>::solveExact() noexcept {
   
   const Type Ns = ( _tf -_t0 )/ _dt ;
    

   std::ofstream fn( filename , std::ios::out);
   
   if(!fn)
   {
      std::cerr << "Error opening file " << filename << " in rhsOdeProblem::solveExact()" << std::endl;      
      exit(-1);
   }
   else
   {
       //std::cout << Ns << std::endl ;   
       Type time = _t0 , yt = _u0 ;
   
       fn << time << ' ' << yt << std::endl ;   
       for(std::size_t i=0 ; i < Ns ; i++)
       {
          time += _dt ; 
          yt = analiticalFunction(time,yt);
      
          fn << time << ' ' << yt << std::endl ;   
       }   
      fn.close();
    }
      
}

  }//ode 
 }//numeric
}//mg 
# endif
