MODULE compute  

      use declarations

      implicit none

      CONTAINS

      function first_find(array, VAL) result (loc)
              !Determines the location of the element in the array with the value given in the VAL argument
              ! return -1 if there is no match
                
              implicit none

              integer, intent(in)             :: array(:)
              integer, intent(in)             :: VAL
              integer                         :: loc

              integer                         :: i, num_elements

              num_elements = size(array)
              loc = -1
              do i = 1, num_elements
                      if (array(i) .EQ. VAL) then
                              loc = i
                              exit
                      end if
              end do

              end function


      function initial(N) result (array_S)
      ! initialize spins randomly 


              integer, intent(in)             :: N     ! number of spins
              integer, dimension(N)           :: array_S
              
              integer                         :: i, loc
              integer, dimension(33)           :: seed
              real                            :: num

              ! get system clock count
              call system_clock(i)
              ! set random number seed
              seed = i
              ! seed the random number generator
              call random_seed(PUT=seed)
          
              do loc = 1, N
                      call random_number(num)
                      if (num < 0.5) then
                              array_S(loc) = -1
                      else
                              array_S(loc) = 1
                      end if
              end do

      end function

      pure function totmag(array_S) result (Mag)
              ! sums the numbers of the input array

              integer, intent(in)             :: array_S(:)
              integer                         :: Mag

              Mag = SUM(array_S)

      end function

      pure function tot_energy(array_S, array_nbr, num_nbr) result (E)
              ! Find energy of an ising configuration

              integer, intent(in)               :: array_S(:)
              integer, intent(in)               :: array_nbr(:,:)
              integer, intent(in)               :: num_nbr(:)
              
              integer                           :: E

              integer                           :: N, i, j, nbr

              E = 0
              N = size(array_S)
              do i = 1, N
                do j = 1, num_nbr(i)
                        nbr = array_nbr(i, j)
                        E = E - array_S(i)*array_S(nbr)
                end do
              end do
              E = E/2

        end function



      subroutine gen_lat(irregular, N, array_sites, array_nbr, num_nbr, numspin, numchain)
              ! Generates a 2D lattice made of 1-D chains. Chains are made of (x,y)
              ! sites with given nn distance.  
              ! Create nearest neighbor table

              logical, intent(in)                      :: irregular
             integer, intent(in)                      :: N
             integer, intent(out)                     :: numspin, numchain

             real(DP), intent(out)                    :: array_sites(:,:) 
             integer, intent(out)                     :: array_nbr(:,:), num_nbr(:)

             !integer, allocatable                     :: border(:)

             real(DP)                                 :: x0, y0, theta, x, y, xp, yp, xn, yn
             integer                                  :: num_p, num_n, i, j, m, nnn, site
             integer, dimension(33)                   :: seed
             real(DP), dimension(dimen)               :: dr
             real(DP)                                 :: dist, num

             !allocate(border(N2))

             array_sites = 0.0
             array_nbr = 0.0
             num_nbr = 0.0
             !border = 0.0
             numspin = 0       ! accumulate till rho is met
             numchain = 0       ! count the number of individual spin chains
             m = 1
             ! Create the spins
              ! get system clock count
              call system_clock(i)
              ! set random number seed
              seed = i
              ! seed the random number generator
              call random_seed(PUT=seed)

             do while (numspin < N)   ! one new chain starts as long as this holds
                      ! generate a random point inside the box
                      ! and a random angle for the chain
                      call random_number(num)
                        if (irregular .EQV. .TRUE.) then 
                                x0 = L*num
                        else
                                x0 = 0.0
                        end if
                      call random_number(num)
                        if (irregular .EQV. .TRUE.) then 
                                y0 = L*num
                        else
                                y0 = (m-1)*a
                        end if
                      call random_number(num)
                        if (irregular .EQV. .TRUE.) then 
                                theta = PI*num
                        else 
                                theta = 0.0
                        end if
                        m = m + 1
                      ! increase chain size by 1
                      numchain = numchain + 1
                      ! position of the current spin
                      x = x0
                      y = y0
                      i = 0
                      ! index the final point in the negative direction
                      ! is the same as the index of (x0, y0), in case the next -ve point is out of bounds
                      num_n = numspin + 1       ! record the spin index of (x0, y0)
                      ! go forward on the chain till you hit the edge
                      do while ((x >= 0.0) .AND. (x < L) .AND. (y >= 0.0) .AND. (y < L))
                              ! increase number of spins by one
                              numspin = numspin + 1
                              array_sites(numspin,:) = (/x,y/)
                              ! spin on the border is the last current spin
                              xp = x
                              yp = y
                              ! generate next spin
                              x = x + a*cos(theta)
                              y = y + a*sin(theta)
                      end do
                      ! save the index of the final point
                      num_p = numspin  
                      ! go backwards on the chain
                      x = x0 - a*cos(theta)
                      y = y0 - a*sin(theta)
                      !i = 1
                      do while ((x >= 0.0) .AND. (x < L) .AND. (y >= 0.0) .AND. (y < L))
                              ! increase number of spins by one
                              numspin = numspin + 1
                              array_sites(numspin,:) = (/x,y/)
                              ! spin on the border is the last current spin
                              xn = x
                              yn = y
                              !i = i + 1
                              ! generate next spin
                              x = x - a*cos(theta)
                              y = y - a*sin(theta)
                      end do
                      if (num_p .NE. numspin) then  ! this means the previous while loop was enetered 
                              num_n = numspin    ! The final point in the negative direction
                      end if  
                      ! fill the nearest neighbor array with endpoint information
                      array_nbr(num_p,1) = num_n
                      array_nbr(num_n,1) = num_p
                      !border(num_p) = 1
                      !border(numspin) = 1
               
              end do

              ! Now fill the rest of the nearest neighbor array

              do site = 1, numspin
                      if (array_nbr(site,1) .GT. 0) then
                              nnn = 1  ! number of nearest neighbor of site 'site'
                      else
                              nnn = 0
                      end if
                      do j = 1, numspin
                              dr = array_sites(site,:) - array_sites(j,:)
                              dist = sqrt(dr(1)**2 + dr(2)**2)
                              if ((site .NE. j) .AND. (dist .LE. (a + 0.0001))) then
                                      nnn = nnn + 1
                                      array_nbr(site,nnn) = j
                              end if
                      end do
                      num_nbr(site) = nnn
              end do

      end subroutine

      subroutine wolff(nsteps, array_S, array_nbr, num_nbr, p, binder, mean_E, mean_E2)

              ! run wolff algorithm. Return binder ratio

              integer, intent(in)             :: nsteps
              ! Spin array, array of neighbors, array of neighbor numbers
              integer, intent(inout)          :: array_S(:)
              integer, intent(in)             :: array_nbr(:,:), num_nbr(:)
              ! "magic" value for 100% acceptance for flipping domains
              real(DP), intent(in)            :: p
              real(DP), intent(out)           :: binder, mean_E, mean_E2

              integer, allocatable            :: pocket(:), cluster(:)
              ! total magnetization for each step
              real(DP), dimension(nsteps)     :: M, E, E2

              integer                         :: i, j, k, l, step, nn, N
              integer, dimension(33)           :: seed
              real(DP)                        :: num, mean_m, mean_m2, mean_m4, R

              N = size(array_S)
              allocate(pocket(N), cluster(N))
              ! set the pocket array to zero. This means pocket is empty.  
              pocket = 0
              ! set the cluster array to 1. This means the cluster is empty. 
              cluster = 1
              M = 0

              ! get system clock count
              call system_clock(i)
              ! set random number seed
              seed = i
              ! seed the random number generator
              call random_seed(PUT=seed)
              
              do step = 1, nsteps
                      ! generate random number [0,1)
                      call random_number(num)
                      ! map from [0,1) to {1,2,3,...N}
                      ! j = m + FLOOR((n+1-m)*u)  ! We want to choose one from n-m+1 integers
                      k = 1 + FLOOR(N*num)
                      ! put site k in the pocket by turning the k-th index from 0 to 1. 
                      pocket(k) = 1
                      ! put site k in the cluster by turning the k-th index from 1 to -1. 
                      cluster(k) = -1
                      ! run till any of the pocket site is occupied
                      do while (first_find(pocket, 1) >=0)
                              ! the first occurance of value in an array
                              j = first_find(pocket, 1)
                              ! take the site out of the pocket by turning the pocket array value to zero.
                              pocket(j) = 0
                              ! loop over neighbors of the site you just took out of the pocket
                              do l = 1, num_nbr(j)
                                      nn = array_nbr(j,l)
                                      ! if the spins are aligned AND if the nn is not in the cluster already
                                      if ((array_S(nn) .EQ. array_S(j)) .AND. (cluster(nn) .EQ. 1)) then
                                              call random_number(num)
                                              ! if the random number drawn is less than p
                                              if (num .LE. p) then
                                                      ! include the nn in the pocket
                                                      pocket(nn) = 1
                                                      ! include the nn in the cluster
                                                      cluster(nn) = -1
                                              end if
                                      end if
                              end do
                      end do
                      ! at this stage a cluster is built and the pocket is empty
                      array_S = array_S*cluster
                      cluster = 1
                      M(step) = totmag(array_S)
                      E(step) = tot_energy(array_S, array_nbr, num_nbr)
                      E2(step) = (E(step))**2

              end do
              mean_E = sum(E)/(N*nsteps)
              mean_E2 = sum(E2)/(N*nsteps)
              mean_m = sum(M)/(N*nsteps)
              mean_m4 = sum(M**4)/nsteps
              mean_m2 = sum(M**2)/nsteps
              R = mean_m4/(mean_m2**2)
              binder = 0.5*(3-R)

              deallocate(pocket, cluster)
      end subroutine


END MODULE compute
