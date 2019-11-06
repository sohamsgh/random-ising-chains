      module declarations

      ! Store constants of the system being calculated

      implicit none

      save

      ! DP
      integer DP
      parameter (DP=kind(1.d0))

      ! PI
      real(DP) PI
      parameter (PI=4.D0*DATAN(1.D0))
      
      ! Boltzmann's constant
      real(DP) Kb
      parameter (Kb=1.0)

      ! Dimension
      integer dimen
      parameter (dimen=2)
      
      ! distance between points on a chain 
      real(DP) a
      parameter (a=1.0)

      
      !linear dimension of box.
      integer L
      parameter (L=8)

      ! density of spins
      real(DP) rho
      parameter (rho=1.0)

      ! circular or square spins
      logical circular  
      parameter (circular=.FALSE.)

      ! Temperature beginning, end
      real(DP) T1, T2
      parameter (T1=1.0,T2=4.0)
      
      ! Temperature steps.
      integer ntemp
      parameter(ntemp=30)

      ! number of thermalization steps
      integer nthermal
      parameter (nthermal=100000)
      
      ! number of steps
      integer nsteps
      parameter (nsteps=1000000)
      
      logical irregular
      parameter (irregular=.FALSE.)
      
      end module declarations

