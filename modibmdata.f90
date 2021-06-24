!> \file modibmdata.f90
!! By Michael Koene (MK), TU Delft, section Atmospheric Physics, 28 January 2019
!! cstep: subroutines airslabsum and wallslabsum moved from modmpi to here to avoid mutual dependencies

module modibmdata

  implicit none
  save
  public :: applydamping,airslabsum,wallslabsum

  logical :: lapply_ibm     = .false.        !< Switch to enable immersed boundary method 
  logical :: lreadfile_obstacles  = .false.  !< Switch to read positions and model level height from a file
  logical :: lnoslip        = .false.        !< Switch to use a no slip condition for the walls
  logical :: lwallfunc      = .true.         !< Switch to use wallfunctions to describe wallshear
  logical :: lfluxform      = .true.         !< Switch to use the fluxform of the advection in advection correction
  logical :: lpoislast      = .true.         !< Switch to use the Poisson solver after the Immersed boundary method
                                             !  .false. will set the order to: ZeroVelocity -> PoissonSolver -> IBM

  real    :: thlwall        = 293.           !< Wall temperature for temperature flux at the sides and top of the buildings
  real    :: ct             = 1.             !< coefficient for temperature flux near the wall 

  !< Field for the wall shear stress
  real, allocatable    :: damping(:,:,:)
  !< Field for the immersed boundary height
  integer, allocatable :: bc_height(:,:)     !< Height of immersed boundary at grid pos x,y
  !< Boolean for applying IBM
  logical, allocatable :: libm (:,:,:)       
  !< Number of grid points in a slab excluding obstacles, and the number of obstacle points
  integer, allocatable :: Nair (:), Nwall(:)
 

contains

  !MK: Subroutine airslabsum - Calculates the slabsum of the cells excluding the immersed boundary cells
  !                            Nair is the number of air cells at each height usable for calculating correct averages.
  !                            it is used in modthermodynamics and modgenstat
  subroutine airslabsum(aver,ks,kf,var,ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes) !cstep ,Nair)
    use modglobal, only: i1,j1
    use modfields, only: ksfc !cstep IBM
    use modmpi, only: my_real,mpi_sum,comm3d,mpierr !cstep IBM ,myidx,myidy
    implicit none

    !cstep integer, intent(out) :: Nair(ks:kf)
    integer :: ks,kf
    integer :: ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes
    real    :: aver(ks:kf)
    real    :: var (ib:ie,jb:je,kb:ke)
    real    :: averl(ks:kf)
   !cstep real    :: avers(ks:kf)
    integer :: i,j,k !cstep IBM ,kmin

    averl(:)    = 0.
    aver(:)    = 0.

    do i=2,i1
    do j=2,j1
    !cstep IBM kmin=bc_height(i+myidx*imax,j+myidy*jmax)+1
    !cstep IBM  do k=kmin,kes
    do k=ksfc(i,j),kes
      averl(k) = averl(k)+var(i,j,k)
    enddo
    enddo
    enddo

    call MPI_ALLREDUCE(averl, aver, kf-ks+1,  MY_REAL, &
                       MPI_SUM, comm3d,mpierr)
    aver = aver/Nair

    return
  end subroutine airslabsum

  !MK: Subroutine wallslabsum - Calculates the slabsum of the cells containing buildings ignoring the air cells
  !                            Nwall is the number of immersed boundary cells in total usable for calculating correct averages.
  subroutine wallslabsum(aver,ks,kf,var,ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes) !cstep ,Nwall)
    use modfields, only: ksfc !cstep IBM
    use modglobal, only: i1,j1,k1 !imax, jmax
    use modmpi, only: my_real,mpi_integer,mpi_sum,comm3d,mpierr !cstep IBM ,myidx,myidy 

    implicit none

    !cstep integer, intent(out) :: Nwall
    integer :: ks,kf
    integer :: ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes
    real    :: aver(kb:ke)
    real    :: var (ib:ie,jb:je,kb:ke)
!cstep    real    :: averl
!cstep    real    :: avers
!cstep    integer :: Nwalll , computation done only once in modibm
!cstep    integer :: Nwalls
    integer :: i,j,k  !kmin  !cstep IBM
    real,allocatable, dimension(:):: averl

    allocate( averl(k1))

    averl(:)    = 0.
!cstep    avers    = 0.
!cstep    aver(:)     = 0.
!cstep    Nwalll    = 0
!cstep    Nwalls    = 0
!cstep    Nwall     = 0

    do i=2,i1
    do j=2,j1
    !cstep IBM kmin=bc_height(i+myidx*imax,j+myidy*jmax)+1
    !cstep IBM do k=1,kmin-1
    do k=1,ksfc(i,j) -1
      averl(k) = averl(k)+var(i,j,k)
 !cstep     Nwalll = Nwalll+1
    enddo
    enddo
    enddo

!cstep    call MPI_ALLREDUCE(averl, avers, 1,  MY_REAL, &
!cstep                       MPI_SUM, comm3d,mpierr)
     call MPI_ALLREDUCE(averl, aver, k1,  MY_REAL, MPI_SUM, comm3d,mpierr)
!cstep    call MPI_ALLREDUCE(Nwalll, Nwalls, 1,  MPI_INTEGER, &
!cstep                       MPI_SUM, comm3d,mpierr)

    !cstep aver = aver + avers
    !cstep Nwall = Nwall + Nwalls

    deallocate(averl)

    return
  end subroutine wallslabsum


 
  subroutine applydamping(ekm,i,j,k)
    
    implicit none
    real, intent(inout) :: ekm
    integer, intent(in) :: i,j,k

    if ((.not. (lapply_ibm)).or.(.not.(lwallfunc))) return
    ekm=ekm*damping(i-1,j-1,k)
    return
  end subroutine applydamping
end module modibmdata
