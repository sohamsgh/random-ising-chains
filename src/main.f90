program main

!################################################################################
!
!       Compute the binder ratio for a 2-D system of random Ising chains 
!
!                               SSG, ZL, RC
!                                 Oct 2019      
!
!################################################################################
        use declarations
        use compute
        implicit none
  
        integer                                   :: N, N2, numspin,numchain
  
  
        integer, allocatable                      :: array_spins(:) 
        real(DP), allocatable                     :: array_sites(:,:)
  
        real(DP), allocatable                     :: tmp_array_sites(:,:)
        integer, allocatable                      :: array_nbr(:,:), num_nbr(:)
        integer, allocatable                      :: tmp_array_nbr(:,:),tmp_num_nbr(:)
  
        real(DP)                                  :: binder
        integer                                   :: i
        character(len=256)                      :: outfile 
        real(DP)                                  :: p, T, beta, mean_E, mean_E2, Cv

  
        ! Number of spins
        if  (circular .EQV. .TRUE.) then
                N = int(rho*(L**2)*(4.0/PI))
        else
                N = int(rho*(L**2))
        end if
        ! Since a chain than has begun is allowed to be constructed till it is complete
        ! there can be at the most N + L*sqrt(2) +2 spins after N is defined
        ! This is when the last chain is along the diagonal
        ! Extra two for spins on the edges
        N2 = N + int(L*sqrt(2.0)) +2
        !the actual coordinates of the sites
        allocate(tmp_array_sites(N2, dimen))
        ! the nearest neighbor array
        ! a spin can have at the most N2 neighbors
        allocate(tmp_array_nbr(N2, N2), tmp_num_nbr(N2))
  
        ! initialize
        tmp_array_sites = 0.0
        tmp_array_nbr = 0
        tmp_num_nbr = 0
        
        call getarg(1, outfile)

        call gen_lat(irregular, N, tmp_array_sites, tmp_array_nbr, tmp_num_nbr, numspin, numchain)
  
        ! after getting numspin, allocate the arrays
        allocate(array_sites(numspin, dimen))
        allocate(array_nbr(numspin, numspin), num_nbr(numspin))
        allocate(array_spins(numspin))

        ! initialize
        array_sites = 0.0
        array_nbr = 0
        num_nbr = 0
        array_spins = 0
  
        ! create the arrays
        array_sites = tmp_array_sites(1:numspin, :)
        array_nbr = tmp_array_nbr(1:numspin, 1:numspin)
        num_nbr = tmp_num_nbr(1:numspin)
        ! create initial spins
        array_spins = initial(numspin)
  
        ! deallocate temporary arrays
        deallocate(tmp_array_sites, tmp_array_nbr, tmp_num_nbr)
        
        ! Write site information
        ! write the sites
        open(121, file = 'chains.xyz', status = 'replace')
        write(121, '(A)') "2-D random Ising chains"
        do i = 1, numspin
                write(121, '(2F10.5)') array_sites(i,1), array_sites(i,2)
        end do
        close(121)
        ! Write spins and nearest neighbor information
        do i = 1, numspin
                write(6,'(I5, I5, 10I5)') array_spins(i), num_nbr(i), array_nbr(i,1:10)
        end do

        ! open output file
        open(221, file=trim(adjustl(outfile)), status = 'replace') 
        write(221,'(A, 4X, I6, 4X, F5.3, 4X, F5.3, 4X, I5 )') "#", L, rho, a, numspin

        do i = 1, ntemp
                ! create initial spins
                array_spins = initial(numspin)
                T = T1 + (T2-T1)/float(ntemp)*(i-1)
                beta = 1/(Kb*T)
                p = 1.0 - exp(-2.0*beta)
                call wolff(nthermal, array_spins, array_nbr, num_nbr, p, binder, mean_E, mean_E2)
                call wolff(nsteps, array_spins, array_nbr, num_nbr, p, binder, mean_E, mean_E2)
                Cv = (mean_E2 - (mean_E)**2)*((beta**2)/numspin)
                write(221,'(F5.3, 4x, F10.6)')  T, mean_E
                call FLUSH(221)
        end do
                
        close(221) 
        deallocate(array_sites, array_spins, array_nbr, num_nbr)

  end program main
            
