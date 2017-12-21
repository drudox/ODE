# ifndef __RUNGEKUTTA_H__
# define __RUNGEKUTTA_H__ 

# include "../rhsOdeProblem.H"
# include "../OdeSolver.H"


namespace mg { 
                namespace numeric {
                                    namespace ode {


/*-------------------------------------------------------------------------------
 *    
 *    Base class for all the Runge-Kutta schemes 
 *    dy/dt = f(y,t)
 *
 *    @Marco Ghiani October 2017, Glasgow UK
 *
 ------------------------------------------------------------------------------*/


template <typename Type = double> 
class RungeKutta :      
                        public OdeSolver<Type>
{
    
   public: 
      
      RungeKutta(const rhsOdeProblem<Type>& that ) noexcept : 
                                                                OdeSolver<Type>{that} 
      {
         Nstep = (that.tf()-that.t0())/that.dt();                 
         reSize();   
      }                                
      
      virtual ~RungeKutta() = default ;
      
      using OdeSolver<Type>::rhs;
      
    //  virtual void solve(const std::string& ) override = 0 ;
    //  virtual void solve() noexcept override           = 0 ;
      virtual void solve(const std::string filename) override = 0;
      virtual void solve() noexcept  override                 = 0;
       
      void reSize() noexcept ; 

   protected:
      
      std::size_t Nstep ; 
      
      std::vector<Type> uc ;

      Type k1;
      Type k2;
      Type k3;
      Type k4;
      Type K5;
      Type k6; 
      Type k7;
      Type k8; // support RK 8-th order accuracy schemes 

};

template<typename T>
void RungeKutta<T>::reSize() noexcept 
{
    uc.resize(Nstep+1);
}




  }//ode 
 }//numeric
}//mg
# endif 
