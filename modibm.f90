!> \file modibm.f90
!! By Michael Koene, email: mich4all@live.nl, TU Delft, section Atmospheric Physics, October 8, 2019
!! TODO: Write comments to output file
!!       clean up code (write statements)
!!       Test performance

module modibm

  use modibmdata, only : lapply_ibm,lpoislast, lreadfile_obstacles, lnoslip, lwallfunc, lfluxform, damping, thlwall, ct,&
                           bc_height,libm,Nair,Nwall
  implicit none
  save
  public :: initibm, exitibm,&
            applyibm, zerowallvelocity

  ! Fields
!libm  logical, allocatable :: limmersed_boundary(:,:,:)        !< Boolean array where .true. is the immersed boundary

  !< Shear layers in x,y,z-directions
  logical, allocatable :: lshear_x(:,:,:), lshear_y(:,:,:), lshear_z(:,:,:)
  !< Normal immersed boundary layers for incoming x,y,z-velocities
  logical, allocatable :: lnorm_x(:,:,:), lnorm_y(:,:,:), lnorm_z(:,:,:)

  !< Field for the minimal distance between a cellcenter and a wall
  real, allocatable    :: mindist(:,:,:)

  !< Field for the wall shear stress
  real, allocatable    :: shear(:,:,:,:)

  !< Fields for velocities at the ghost positions (k=0)
  real, allocatable    :: u0g(:,:),  v0g(:,:)   !cstep, umg(:,:), vmg(:,:)  were computed but are not used

  !< Field for the ekm at the ghost positions (k=0)
  real, allocatable    :: ekmg(:,:)

  real, allocatable    :: tempsvp(:,:,:,:)
  real, allocatable    :: tempthlp(:,:,:)
  real, allocatable    :: tempup(:,:,:), tempvp(:,:,:), tempwp(:,:,:)

!cstep  integer, allocatable :: Nair(:)
!cstep  integer, allocatable :: bc_height(:,:)

contains
  subroutine initibm
    use modglobal,  only : itot, jtot, ih, i1, jh, j1, k1, imax, jmax, kmax, cexpnr, ifnamopt, ifinput, &
                           fname_options, nsv, e12min, cu, cv, ijtot
    use modmpi,     only : myid, mpi_logical, mpi_sum,comm3d, mpierr, MPI_INTEGER, myidx, myidy, my_real, boolexcjs
    use modsurfdata, only : thls
    use modfields,   only : thl0,thlm,qt0,qtm,u0,um,v0,vm,w0,wm,e120,e12m,ksfc
    use modsubgrid,  only : ekh,ekm
    implicit none

    integer, allocatable :: Nairl(:)
    integer       :: i, j, k, ierr  !cstep , kmin
    character(600) :: readstring

    namelist/NAMIBM/ lapply_ibm, lreadfile_obstacles, &
                               lwallfunc, lnoslip, lfluxform, &
                               thlwall, lpoislast, ct


    if(myid==0) then    !first myid
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMIBM,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMIBM'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMIBM'
      endif
      write(6 ,NAMIBM)
      close(ifnamopt)
    endif

    !if(.not.(myid==0)) return

    call MPI_BCAST(lapply_ibm   ,    1, mpi_logical , 0, comm3d, mpierr)

    if (.not. (lapply_ibm)) return


    call MPI_BCAST(lwallfunc        ,    1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(lnoslip          ,    1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(lreadfile_obstacles        ,    1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(lfluxform        ,    1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(thlwall          ,    1, my_real     , 0, comm3d, mpierr)
    call MPI_BCAST(lpoislast        ,    1, mpi_logical , 0, comm3d, mpierr)
    call MPI_BCAST(ct               ,    1, my_real     , 0, comm3d, mpierr)

    if (lnoslip .and. lwallfunc) then
      if(myid==0) print *, 'Problem in namoptions NAMIBM'
      if(myid==0) print *, 'Cannot use both no slip conditions and wall functions for the shear'
      if(myid==0) print *, 'Either set lnoslip to true or lwallfunc to true but not both.'
      if(myid==0) print *, 'lwallfunc = .true. is recommended'
      stop 'ERROR: Problem in namoptions NAMIBM'
    endif

    if (.not. (lnoslip .or. lwallfunc)) then
      if(myid==0) print *, 'Problem in namoptions NAMIBM'
      if(myid==0) print *, 'Cannot go without no slip conditions or wall functions for the shear'
      if(myid==0) print *, 'Either set lnoslip to true or lwallfunc to true by the use of lnoslip = .true. or lwallfunc = .true.'
      if(myid==0) print *, 'lwallfunc = .true. is recommended'
      stop 'ERROR: Problem in namoptions NAMIBM'
    endif

    if (abs(cu)>1e-15 .or. abs(cv)>1e-15) then
      if(myid==0) print *, 'Problem in namoptions'
      if(myid==0) print *, 'cu or cv cannot be nonzero while using IBM'
      if(myid==0) print *, 'The buildings would move in that case'
      if(myid==0) print *, 'Set cu and cv to 0. to solve this problem or simulate without buildings'
      stop 'ERROR: Problem in namoptions NAMIBM with cu and cv'
    endif

    write(6,*) 'allocating fields in modibm'

!libm    allocate(limmersed_boundary (itot+1,jtot+1,k1))
    allocate(bc_height (itot+1,jtot+1))
    allocate(libm(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(Nair(k1))
    allocate(Nairl(k1))
    allocate(Nwall(k1))
    if (lnoslip) then
      !libm allocate(lshear_x (2-ih:imax+ih,2-jh:jmax+jh,kmax))
      !libm allocate(lshear_y (2-ih:imax+ih,2-jh:jmax+jh,kmax))
      !libm allocate(lshear_z (2-ih:imax+ih,2-jh:jmax+jh,kmax))
      allocate(lshear_x (2-ih:i1+ih,2-jh:j1+jh,k1)) !cstep replace by (1:i2,1:j2,k1) (one neighboorpoint ih=jh=1, and adapt boolexcjs)
      allocate(lshear_y (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(lshear_z (2-ih:i1+ih,2-jh:j1+jh,k1))
    endif
    !libm allocate(lnorm_x (1:i1,1:j1,kmax))
    !libm allocate(lnorm_y (1:i1,1:j1,kmax))
    !libm allocate(lnorm_z (1:i1,1:j1,kmax))
    allocate(lnorm_x (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(lnorm_y (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(lnorm_z (2-ih:i1+ih,2-jh:j1+jh,k1))

    if (lwallfunc) then
      allocate(shear(2-ih:i1+ih,2-jh:j1+jh,k1,12))
      allocate(damping(2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(u0g(2-ih:i1+ih,2-jh:j1+jh))
!cstep      allocate(umg(2-ih:i1+ih,2-jh:j1+jh))
      allocate(v0g(2-ih:i1+ih,2-jh:j1+jh))
!cstep      allocate(vmg(2-ih:i1+ih,2-jh:j1+jh))
      allocate(ekmg(2-ih:i1+ih,2-jh:j1+jh))
      allocate(tempsvp(2-ih:i1+ih,2-jh:j1+jh,k1,nsv))
      allocate(tempthlp(2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(tempup(2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(tempvp(2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(tempwp(2-ih:i1+ih,2-jh:j1+jh,k1))
    endif

    write(6,*) 'succesfully allocated fields in modibm'

!libm    limmersed_boundary(:,:,:) = .false.
    bc_height(:,:) = 0
    libm (:,:,:) = .false.
    
    if (lwallfunc) then
      damping(:,:,:)=1.
      shear(:,:,:,:)=0.
      u0g(:,:)=0.
!cstep      umg(:,:)=0.
      v0g(:,:)=0.
!cstep      vmg(:,:)=0.
    endif


    ! Definition of obstacles
    if (myid==0) then
      if (lreadfile_obstacles) then  !< Profile prescribed by use in the file ibm.inp.<expnr>
        write(6,*) 'Reading inputfile in modibm'
        open (ifinput,file='ibm.inp.'//cexpnr)
          do i=1,itot
            read (ifinput,'(a530)') readstring
            !< If the program is unable to read the full line of points increasing the length of the string (a400) might help

            do while (readstring(1:1)=='#')  ! Skip the lines that are commented (like headers)
              read (ifinput,'(a530)') readstring
            end do
            read(readstring,*) (bc_height(i+1,j+1),j=1,jtot)
          end do
        close(ifinput)

        bc_height(1,:)=bc_height(itot+1,:)
        bc_height(:,1)=bc_height(:,jtot+1)

!libm        do i=1,itot
!libm          do j=1,jtot
!libm            do k=1,kmax
!libm              if (k.LE.bc_height(i,j)) then
!libm                limmersed_boundary(i,j,k)=.true.
!libm              endif
!libm            end do
!libm          end do
!libm        end do

        write(6,*) 'Succesfully read inputfile in modibm'
      else           !< Simple block in the middle of the domain
        write(6,*) 'Generating standard boundary in modibm'
        bc_height(NINT(itot*0.5):(NINT(itot*0.5)+1),NINT(jtot*0.5):(NINT(jtot*0.5)+1))=NINT(kmax*0.5) 

        bc_height(1,:)=bc_height(itot+1,:)
        bc_height(:,1)=bc_height(:,jtot+1)

!libm        do i=1,itot
!libm          do j=1,jtot
!libm            do k=1,kmax
!libm              if (k.LE.bc_height(i,j)) then
!libm                limmersed_boundary(i,j,k)=.true.
!libm              endif
!libm            end do
!libm          end do
!libm        end do

        write(6,*) 'Succesfully generated immersed boundary in modibm'
      endif
    endif  !myid==0

!libm    if(myid==0) write(6,*) 'itot,jtot,imax,jmax,i1,j1=',itot,jtot,imax,jmax,i1,j1


!libm    limmersed_boundary(itot+1,:,:)=limmersed_boundary(1,:,:)
!libm    limmersed_boundary(:,jtot+1,:)=limmersed_boundary(:,1,:)
    call MPI_BCAST(bc_height,(itot+1)*(jtot+1),MPI_INTEGER ,0,comm3d,mpierr)
    !cstep call MPI_BCAST(limmersed_boundary,itot*jtot*kmax,MPI_LOGICAL ,0,comm3d,mpierr)
!libm    call MPI_BCAST(limmersed_boundary,(itot+1)*(jtot+1)*k1,MPI_LOGICAL ,0,comm3d,mpierr)

    do i=2,i1
      do j=2,j1
        do k=1,kmax
          if (k.LE.bc_height(i+myidx*imax,j+myidy*jmax)) then
            libm (i,j,k) = .true.
          endif
        end do
      end do
    end do

    call boolexcjs( libm  , 2,i1,2,j1,1,k1,ih,jh)



    Nair(:) = 0.
    Nairl(:) = 0
    do i=2,i1
      do j=2,j1
!cstep         kmin=bc_height(i+myidx*imax,j+myidy*jmax)+1
           ksfc(i,j) = bc_height(i+myidx*imax,j+myidy*jmax) + 1 
!libm  write (6,*) myid,i,j,i+myidx*imax,j+myidy*jmax,ksfc(i,j),limmersed_boundary(i+myidx*imax,j+myidy*jmax,1)
!cstep        do k=kmin,k1
        do k=ksfc(i,j),k1 
          Nairl(k) = Nairl(k)+1
        enddo
!libm        do k=1,k1
!libm          libm(i,j,k) = limmersed_boundary(i+myidx*imax,j+myidy*jmax,k)
!libm        enddo
      enddo
    enddo
    write (6,*) 'Nairl ',Nairl
    call MPI_ALLREDUCE(Nairl, Nair, k1,  MPI_INTEGER, MPI_SUM, comm3d,mpierr)
    write (6,*) 'Nair',Nair
    Nwall(:) = ijtot - Nair(:)

    call constructboundarytypes

    write(6,* ) 'start deallocate'
    deallocate (Nairl)
    write(6,*) 'exit initibm'

    !if (lwallfunc) call mindistance
    return
  end subroutine initibm

  subroutine constructboundarytypes   !< Calculate the positions of the different boundary layers in multiple directions

    use modglobal,  only : imax, i1, ih, jmax, j1, jh, kmax, k1
    use modmpi,     only : myid,boolexcjs
    implicit none
    integer i,j,k,ipos,jpos


    if(myid==0) write(6,*) 'Starting constructboundarytypes in modibm, myid =',myid

    !< Fill the layer types with .false.
    if(lnoslip) then
      lshear_x(:,:,:)=.false.
      lshear_y(:,:,:)=.false.
      lshear_z(:,:,:)=.false.
    endif
    lnorm_x(:,:,:)=.false.
    lnorm_y(:,:,:)=.false.
    lnorm_z(:,:,:)=.false.

    !< Find Shear layers perpendicular to the x-direction and normal layers in x-direction
    if(lnoslip) then  !< MK: lnoslip = .false. is recommended
      do k=1,kmax
       !libm do j=2,jmax
       !libm   do i=2,imax-1
          !libm  ipos=i+myidx*imax
          !libm  jpos=j+myidy*jmax
          !libm   if (.not. (limmersed_boundary(ipos,jpos,k)==limmersed_boundary(ipos+1,jpos,k))) then
          !libm    lshear_x(i,j,k)=limmersed_boundary(ipos+1,jpos,k)
          !libm    lshear_x(i+1,j,k)=limmersed_boundary(ipos,jpos,k)
         do i=2,i1
           do j=2,j1
             if (.not. (libm(i,j,k)==libm(i+1,j,k))) then
              lshear_x(i  ,j,k)=libm(i+1,j,k)
              lshear_x(i+1,j,k)=libm(i,j,k)
              lnorm_x(i+1 ,j,k)=.true.
            endif
          end do
        end do
      end do
      if(myid==0) write(6,*) 'Succesfully found shear and normal layers in x-direction'

      !< Find Shear layers perpendicular to the y-direction and normal layers in y-direction
      do k=1,kmax
        !libm do i=2,imax
        !libm  do j=2,jmax-1
         !libm   ipos=i+myidx*imax
         !libm   jpos=j+myidy*jmax
         !libm   if (.not. (limmersed_boundary(ipos,jpos,k)==limmersed_boundary(ipos,jpos+1,k))) then
         !libm     lshear_y(i,j,k)=limmersed_boundary(ipos,jpos+1,k)
         !libm     lshear_y(i,j+1,k)=limmersed_boundary(ipos,jpos,k)
          do i=2,i1
           do j=2,j1
            if (.not. (libm(i,j,k)==libm(i,j+1,k))) then
              lshear_y(i,j,k)=libm(i,j+1,k)
              lshear_y(i,j+1,k)=libm(i,j,k)
              lnorm_y(i,j+1,k)=.true.
            endif
          end do
        end do
      end do

      !< Find Shear layers perpendicular to the z-direction and normal layers in z-direction
      !libm do i=2,imax
      !libm  do j=2,jmax
      do i=2,i1
        do j=2,j1
          do k=1,kmax-1
           !libm ipos=i+myidx*imax
           !libm jpos=j+myidy*jmax
           !libm if (.not. (limmersed_boundary(ipos,jpos,k)==limmersed_boundary(ipos,jpos,k+1))) then
           !libm   lshear_z(i,j,k)=limmersed_boundary(ipos,jpos,k+1)
           !libm   lshear_z(i,j,k+1)=limmersed_boundary(ipos,jpos,k)
            if (.not. (libm(i,j,k)==libm(i,j,k+1))) then
              lshear_z(i,j,k)  =libm(i,j,k+1)
              lshear_z(i,j,k+1)=libm(i,j,k)
              lnorm_z (i,j,k+1)=.true.
            endif
          end do
        end do
      end do
      !libm call boolexcjs( lnorm_x  , 2,imax,2,jmax,1,kmax,ih,jh)
      !libm call boolexcjs( lnorm_y  , 2,imax,2,jmax,1,kmax,ih,jh)
      !libm call boolexcjs( lnorm_z  , 2,imax,2,jmax,1,kmax,ih,jh)
      !libm call boolexcjs( lshear_x  , 2,imax,2,jmax,1,kmax,ih,jh)
      !libm call boolexcjs( lshear_y  , 2,imax,2,jmax,1,kmax,ih,jh)
      !libm call boolexcjs( lshear_z  , 2,imax,2,jmax,1,kmax,ih,jh)
 
      call boolexcjs( lshear_x  , 2,i1,2,j1,1,k1,ih,jh)
      call boolexcjs( lshear_y  , 2,i1,2,j1,1,k1,ih,jh)
      call boolexcjs( lshear_z  , 2,i1,2,j1,1,k1,ih,jh)


      if(myid==0) write(6,*) 'Succesfully found shear and normal layers in all directions'
    elseif(lwallfunc) then    !< MK: lwallfunc = .true. is recommended
      !< Find normal layers in x-direction
      do k=1,kmax
        do j=2,j1
          do i=2,i1
         !libm   ipos=i+myidx*imax
         !libm   jpos=j+myidy*jmax
         !libm   if (.not. (limmersed_boundary(ipos,jpos,k)==limmersed_boundary(ipos-1,jpos,k))) then
            if (.not. (libm(i,j,k)==libm(i-1,j,k))) then
              lnorm_x(i,j,k)=.true.
            endif
          end do
        end do
      end do
      if(myid==0) write(6,*) 'Succesfully found normal layers in x-direction'
      !write(6,*) 'jmax = ',jmax

      !< Find normal layers in y-direction
      do k=1,kmax
        do i=2,i1
          do j=2,j1
           !libm ipos=i+myidx*imax
           !libm jpos=j+myidy*jmax
           !libm if (.not. (limmersed_boundary(ipos,jpos,k)==limmersed_boundary(ipos,jpos-1,k))) then
            if (.not. (libm(i,j,k)==libm(i,j-1,k))) then
              lnorm_y(i,j,k)=.true.
            endif
          end do
        end do
      end do

      !< Find normal layers in z-direction
      !libm do i=1,imax
       do i=2,i1
        do j=2,j1
          do k=2,kmax
           !libm ipos=i+myidx*imax
           !libm jpos=j+myidy*jmax
           !libm if (.not. (limmersed_boundary(ipos,jpos,k)==limmersed_boundary(ipos,jpos,k-1))) then
            if (.not. (libm(i,j,k)==libm(i,j,k-1))) then
              lnorm_z(i,j,k)=.true.
            endif
          end do
        end do
      end do


      if(myid==0) write(6,*) 'Succesfully found normal layers in all directions'
    endif

    call boolexcjs( lnorm_x  , 2,i1,2,j1,1,k1,ih,jh)
    call boolexcjs( lnorm_y  , 2,i1,2,j1,1,k1,ih,jh)
    call boolexcjs( lnorm_z  , 2,i1,2,j1,1,k1,ih,jh)

    if(myid==0) write(6,*) 'finished constructboundarytypes'



  end subroutine constructboundarytypes

!  subroutine mindistance !< Determine the distance between the cellcenter and nearest wall for each cell (i,j,k)
!
!    use modglobal, only : itot, jtot, kmax, zh, zf, dx, dy
!    implicit none
!    integer :: i,j,k,iw,jw,kw
!    real    :: dist
!
!    mindist(:,:,:)=10*(itot+jtot+kmax)
!    write(6,*) 'Calculating minimal distances between each point in space and the walls'
!    do i=1,i1
!      do j=1,j1
!        do k=1,kmax
!          mindist(i,j,k)=zh(k)
!          do iw=1,itot+1
!            do jw=1,jtot+1
!              do kw=1,kmax
!                if(lnorm_x(iw,jw,kw)) then
!                  dist=sqrt(((iw-0.5-i)*dx)**2+((jw-j)*dy)**2+(zh(iw)-zh(i))**2)
!                  if (dist < mindist(i,j,k)) then
!                    mindist(i,j,k) = dist
!                  endif
!                elseif(lnorm_y(iw,jw,kw)) then
!                  dist=sqrt(((iw-i)*dx)**2+((jw-0.5-j)*dy)**2+(zh(iw)-zh(i))**2)
!                  if (dist < mindist(i,j,k)) then
!                    mindist(i,j,k) = dist
!                  endif
!                elseif(lnorm_z(iw,jw,kw))then
!                  dist=sqrt(((iw-i)*dx)**2+((jw-j)*dy)**2+(zf(iw)-zh(i))**2)
!                  if (dist < mindist(i,j,k)) then
!                    mindist(i,j,k) = dist
!                  endif
!                endif
!              end do
!            end do
!          end do
!        end do
!      end do
!    end do
!
!    return
!  end subroutine mindistance

  subroutine exitibm
    use modmpi, only : myid
    implicit none

    !deallocate(bc_height,Nair,Nwall)
      if (.not. (lapply_ibm)) return
    if(myid==0) write(6,*) 'Starting with exitibm'
    deallocate(Nair,Nwall) !Andere locatie, en libm is ook verwijderd 1x, Bart Warmerdam
    deallocate(bc_height)
    if(lnoslip) then
      deallocate(lshear_x)
      deallocate(lshear_y)
      deallocate(lshear_z)
    elseif(lwallfunc) then
      !deallocate(mindist)
      !write(6,*) 'deallocating shear'
      deallocate(shear)
      !write(6,*) 'deallocating damping'
      deallocate(damping)
      !write(6,*) 'deallocating ghost points'
      deallocate(u0g)
!cstep      deallocate(umg)
      deallocate(v0g)
!cstep      deallocate(vmg)
      deallocate(ekmg)
      deallocate(tempsvp)
      deallocate(tempthlp)
      deallocate(tempup)
      deallocate(tempvp)
      deallocate(tempwp)
    endif
    deallocate(lnorm_x)
    deallocate(lnorm_y)
    deallocate(lnorm_z)
!libm    deallocate(limmersed_boundary)
    deallocate(libm)
    if(myid==0) write(6,*) 'Finished with exitibm'
    return
  end subroutine exitibm

  subroutine applyibm(simid)   !< apply immersed boundary method
    use modfields,      only : um, vm, wm, u0, v0, w0, up, vp, wp
    use modglobal,      only : rk3step, itot, imax, jmax, jtot, kmax, i1, j1, k1, ih, jh, dt, rdt, timee, dx, dy, dzh, dzf, ekmin, nsv, e12min
    use modfields,      only : rhobf, rhobh, thl0, thlp, sv0, svp, e12p, thlm, e12m
    use modsubgriddata, only : ekm
    use modmicrodata,   only : nu_a
    use modmpi,         only : myid, excjs !libm ,myidx,myidy
    use modnudgeboundary,        only : Nsim
    
    implicit none
    integer  :: i, j, k, m, nc,ipos,jpos
    real     :: rk3coef,rk3coefi
    real     :: emmo, emom, eomm, empo, epmo, epom, emop, eopm, eomp
    real     :: yplus

    integer  :: maxlocx(3)
integer, intent(in) :: simid
    if (.not. (lapply_ibm .and. simid == Nsim)) return

    if (myid==0) then
      !write(6,*) 'Starting with applyibm, t = ',timee
      !write(6,*) 'The timestep dt = ',dt
    endif

    rk3coef = rdt / (4. - dble(rk3step))
    rk3coefi = 1. / rk3coef

    damping(:,:,:)=1.

    tempsvp(:,:,:,:)=0.
    tempthlp(:,:,:)=0.
    tempup(:,:,:)=0.
    tempvp(:,:,:)=0.
    tempwp(:,:,:)=0.

    if(lnoslip) then  !< MK: lnoslip = .false. is recommended
       do i=2,i1  !1+myid*imax,myid*imax+1+imax
          do j=2,j1              !1+myid*jmax,myid*jmax+1+jmax
             do k=2,kmax          !cstep special treatment for k=1
            !< Shear layer first
               if (lshear_x(i,j,k)) then
                  vp(i,j,k)=-vm(i,j,k)*rk3coefi
                  wp(i,j,k)=-wm(i,j,k)*rk3coefi
               endif
               if (lshear_y(i,j,k)) then
                  up(i,j,k)=-um(i,j,k)*rk3coefi
                  wp(i,j,k)=-wm(i,j,k)*rk3coefi
               endif
               if (lshear_z(i,j,k)) then
                  up(i,j,k)=-um(i,j,k)*rk3coefi
                  vp(i,j,k)=-vm(i,j,k)*rk3coefi
               endif
              enddo !k  
           enddo !j
        enddo ! i
    elseif(lwallfunc) then         !< !< MK: lwallfunc = .true. is recommended
      do i=2,i1  !1+myid*imax,myid*imax+1+imax
        do j=2,j1              !1+myid*jmax,myid*jmax+1+jmax
           do k=2,kmax          !cstep special treatment for k=1

            if (lnorm_x(i,j,k)) then     !< Wall in x-direction

              emmo = 0.25  * ( &
                ekm(i,j,k)+ekm(i,j-1,k)+ekm(i-1,j-1,k)+ekm(i-1,j,k)  )

              call wallaw(v0(i-1,j,k),0.5*dx,nu_a,shear(i-1,j,k,1))
              call wallaw(v0(i,j,k),  0.5*dx,nu_a,shear(i,j,k,2))

              tempvp(i-1,j,k) = tempvp(i-1,j,k) - 0.5 * emmo*((v0(i,j,k)-v0(i-1,j,k))/dx) / dx - 0.5 * shear(i-1,j,k,1)/dx
              tempvp(i,j,k)   = tempvp(i,j,k)   + 0.5 * emmo*((v0(i,j,k)-v0(i-1,j,k))/dx) / dx - 0.5 * shear(i,j,k,2)  /dx
              !vp(i-1,j,k) = vp(i-1,j,k) - 0.5 * emmo*((v0(i,j,k)-v0(i-1,j,k))/dx) / dx - 0.5 * shear(i-1,j,k,1)/dx
              !vp(i,j,k)   = vp(i,j,k)   + 0.5 * emmo*((v0(i,j,k)-v0(i-1,j,k))/dx) / dx - 0.5 * shear(i,j,k,2)  /dx

              empo = 0.25  * ( &
                ekm(i,j+1,k)+ekm(i,j,k)+ekm(i-1,j,k)+ekm(i-1,j+1,k)  )

              call wallaw(v0(i-1,j+1,k),0.5*dx,nu_a,shear(i-1,j+1,k,1))
              call wallaw(v0(i,j+1,k),  0.5*dx,nu_a,shear(i,j+1,k,2))

              tempvp(i-1,j+1,k) = tempvp(i-1,j+1,k) - 0.5 * empo*((v0(i,j+1,k)-v0(i-1,j+1,k))/dx) / dx - 0.5 * shear(i-1,j+1,k,1)/dx
              tempvp(i,j+1,k)   = tempvp(i,j+1,k)   + 0.5 * empo*((v0(i,j+1,k)-v0(i-1,j+1,k))/dx) / dx - 0.5 * shear(i,j+1,k,2)  /dx
              !vp(i-1,j+1,k) = vp(i-1,j+1,k) - 0.5 * empo*((v0(i,j+1,k)-v0(i-1,j+1,k))/dx) / dx - 0.5 * shear(i-1,j+1,k,1)/dx
              !vp(i,j+1,k)   = vp(i,j+1,k)   + 0.5 * empo*((v0(i,j+1,k)-v0(i-1,j+1,k))/dx) / dx - 0.5 * shear(i,j+1,k,2)  /dx

              emom = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,k-1) + ekm(i-1,j,k-1) ) ) / &
                    ( 4.   * dzh(k) )

              call wallaw(w0(i-1,j,k),0.5*dx,nu_a,shear(i-1,j,k,3))
              call wallaw(w0(i,j,k),  0.5*dx,nu_a,shear(i,j,k,4))

              tempwp(i-1,j,k) = tempwp(i-1,j,k) - 0.5 * emom * ((w0(i,j,k)-w0(i-1,j,k))/dx)/dx - 0.5 * shear(i-1,j,k,3)/dx
              tempwp(i,j,k)   = tempwp(i,j,k)   + 0.5 * emom * ((w0(i,j,k)-w0(i-1,j,k))/dx)/dx - 0.5 * shear(i,j,k,4)  /dx
              !wp(i-1,j,k) = wp(i-1,j,k) - 0.5 * emom * ((w0(i,j,k)-w0(i-1,j,k))/dx)/dx - 0.5 * shear(i-1,j,k,3)/dx/dx
              !wp(i,j,k)   = wp(i,j,k)   + 0.5 * emom * ((w0(i,j,k)-w0(i-1,j,k))/dx)/dx - 0.5 * shear(i,j,k,4)  /dx/dx

              emop = ( dzf(k) * ( ekm(i,j,k+1)  + ekm(i-1,j,k+1)  )  + &
                      dzf(k+1)  * ( ekm(i,j,k) + ekm(i-1,j,k) ) ) / &
                    ( 4.   * dzh(k+1) )

              call wallaw(w0(i-1,j,k+1),0.5*dx,nu_a,shear(i-1,j,k+1,3))
              call wallaw(w0(i,j,k+1),  0.5*dx,nu_a,shear(i,j,k+1,4))

              tempwp(i-1,j,k+1) = tempwp(i-1,j,k+1) - 0.5 * emop * ((w0(i,j,k+1)-w0(i-1,j,k+1))/dx)/dx - 0.5 * shear(i-1,j,k+1,3)/dx
              tempwp(i,j,k+1)   = tempwp(i,j,k+1)   + 0.5 * emop * ((w0(i,j,k+1)-w0(i-1,j,k+1))/dx)/dx - 0.5 * shear(i,j,k+1,4)  /dx
              !wp(i-1,j,k+1) = wp(i-1,j,k+1) - 0.5 * emop * ((w0(i,j,k+1)-w0(i-1,j,k+1))/dx)/dx - 0.5 * shear(i-1,j,k+1,3)/dx/dx
              !wp(i,j,k+1)   = wp(i,j,k+1)   + 0.5 * emop * ((w0(i,j,k+1)-w0(i-1,j,k+1))/dx)/dx - 0.5 * shear(i,j,k+1,4)  /dx/dx

              call xwallscalar(i,j,k,thl0,tempthlp)
              call xwalltemp(i,j,k,thlm,tempthlp)

              do nc=1,nsv
                call xwallscalar(i,j,k,sv0(:,:,:,nc),tempsvp(:,:,:,nc))
              end do
              call xwalle12(i,j,k)

            endif
            if (lnorm_y(i,j,k)) then     !< Wall in y-direction
              emmo = 0.25  * ( &
                ekm(i,j,k)+ekm(i,j-1,k)+ekm(i-1,j-1,k)+ekm(i-1,j,k)  )

              call wallaw(u0(i,j-1,k),0.5*dy,nu_a,shear(i,j-1,k,5))
              call wallaw(u0(i,j,k)  ,0.5*dy,nu_a,shear(i,j,k,6))

              tempup(i,j-1,k) = tempup(i,j-1,k) - 0.5 * emmo * ((u0(i,j,k)-u0(i,j-1,k))/dy)/dy - 0.5 * shear(i,j-1,k,5)/dy
              tempup(i,j,k)   = tempup(i,j,k)   + 0.5 * emmo * ((u0(i,j,k)-u0(i,j-1,k))/dy)/dy - 0.5 * shear(i,j,k,6)  /dy
              !up(i,j-1,k) = up(i,j-1,k) - 0.5 * emmo * ((u0(i,j,k)-u0(i,j-1,k))/dy)/dy - 0.5 * shear(i,j-1,k,5)/dy
              !up(i,j,k)   = up(i,j,k)   + 0.5 * emmo * ((u0(i,j,k)-u0(i,j-1,k))/dy)/dy - 0.5 * shear(i,j,k,6)  /dy

              epmo = 0.25  * ( &
                ekm(i+1,j,k)+ekm(i+1,j-1,k)+ekm(i,j-1,k)+ekm(i,j,k)  )

              call wallaw(u0(i+1,j-1,k),0.5*dy,nu_a,shear(i+1,j-1,k,5))
              call wallaw(u0(i+1,j,k)  ,0.5*dy,nu_a,shear(i+1,j,k,6))

              tempup(i+1,j-1,k) = tempup(i+1,j-1,k) - 0.5 * epmo * ((u0(i+1,j,k)-u0(i+1,j-1,k))/dy)/dy - 0.5 * shear(i+1,j-1,k,5)/dy
              tempup(i+1,j,k)   = tempup(i+1,j,k)   + 0.5 * epmo * ((u0(i+1,j,k)-u0(i+1,j-1,k))/dy)/dy - 0.5 * shear(i+1,j,k,6)  /dy
              !up(i+1,j-1,k) = up(i+1,j-1,k) - 0.5 * epmo * ((u0(i+1,j,k)-u0(i+1,j-1,k))/dy)/dy - 0.5 * shear(i+1,j-1,k,5)/dy
              !up(i+1,j,k)   = up(i+1,j,k)   + 0.5 * epmo * ((u0(i+1,j,k)-u0(i+1,j-1,k))/dy)/dy - 0.5 * shear(i+1,j,k,6)  /dy

              eomm = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i,j-1,k)  )  + &
                dzf(k) * ( ekm(i,j,k-1) + ekm(i,j-1,k-1) ) ) / ( 4.  * dzh(k) )

              call wallaw(w0(i,j-1,k),0.5*dy,nu_a,shear(i,j-1,k,7))
              call wallaw(w0(i,j,k)  ,0.5*dy,nu_a,shear(i,j,k,8))

              tempwp(i,j-1,k) = tempwp(i,j-1,k) - 0.5 * eomm * ((w0(i,j,k)-w0(i,j-1,k))/dy)/dy - 0.5 * shear(i,j-1,k,7)/dy
              tempwp(i,j,k)   = tempwp(i,j,k)   + 0.5 * eomm * ((w0(i,j,k)-w0(i,j-1,k))/dy)/dy - 0.5 * shear(i,j,k,8)  /dy
              !wp(i,j-1,k) = wp(i,j-1,k) - 0.5 * eomm * ((w0(i,j,k)-w0(i,j-1,k))/dy)/dy - 0.5 * shear(i,j-1,k,7)/dy
              !wp(i,j,k)   = wp(i,j,k)   + 0.5 * eomm * ((w0(i,j,k)-w0(i,j-1,k))/dy)/dy - 0.5 * shear(i,j,k,8)  /dy

              eomp = ( dzf(k) * ( ekm(i,j,k+1)  + ekm(i,j-1,k+1)  )  + &
                dzf(k+1) * ( ekm(i,j,k) + ekm(i,j-1,k) ) ) / ( 4.  * dzh(k+1) )

              call wallaw(w0(i,j-1,k+1),0.5*dy,nu_a,shear(i,j-1,k+1,7))
              call wallaw(w0(i,j,k+1)  ,0.5*dy,nu_a,shear(i,j,k+1,8))

              tempwp(i,j-1,k+1) = tempwp(i,j-1,k+1) - 0.5 * eomp * ((w0(i,j,k+1)-w0(i,j-1,k+1))/dy)/dy - 0.5 * shear(i,j-1,k+1,7)/dy
              tempwp(i,j,k+1)   = tempwp(i,j,k+1)   + 0.5 * eomp * ((w0(i,j,k+1)-w0(i,j-1,k+1))/dy)/dy - 0.5 * shear(i,j,k+1,8)  /dy
              !wp(i,j-1,k+1) = wp(i,j-1,k+1) - 0.5 * eomp * ((w0(i,j,k+1)-w0(i,j-1,k+1))/dy)/dy - 0.5 * shear(i,j-1,k+1,7)/dy
              !wp(i,j,k+1)   = wp(i,j,k+1)   + 0.5 * eomp * ((w0(i,j,k+1)-w0(i,j-1,k+1))/dy)/dy - 0.5 * shear(i,j,k+1,8)  /dy

              call ywallscalar(i,j,k,thl0,tempthlp)
              call ywalltemp(i,j,k,thlm,tempthlp)
              do nc=1,nsv
                call ywallscalar(i,j,k,sv0(:,:,:,nc),tempsvp(:,:,:,nc))
              end do
              call ywalle12(i,j,k)

            endif
            if (lnorm_z(i,j,k)) then     !< Wall in z-direction

            endif
            !TODO: Maybe remove mindist and put in 0.5*dx,0.5*dy,0.5*dz depending on direction
            if (lnorm_x(i,j,k) .or. lnorm_y(i,j,k) .or. lnorm_z(i,j,k)) then !< Calculate and apply Piomelli wall damping
              !yplus= mindist(i,j,k)*sqrt(sum(abs(shear(i,j,k,:))))/nu_a
              if(lnorm_x(i,j,k)) then
                yplus = 0.5 * dx * sqrt(sum(abs(shear(i,j,k,1:4))))/nu_a
                damping(i,j,k)   = min(damping(i,j,k),  1.-exp(-(yplus*0.04)**3.))
                yplus = 0.5 * dx * sqrt(sum(abs(shear(i-1,j,k,1:4))))/nu_a
                damping(i-1,j,k) = min(damping(i-1,j,k),1.-exp(-(yplus*0.04)**3.))
              endif
              if(lnorm_y(i,j,k)) then
                yplus = 0.5 * dy * sqrt(sum(abs(shear(i,j,k,5:8))))/nu_a
                damping(i,j,k)   = min(damping(i,j,k),  1.-exp(-(yplus*0.04)**3.))
                yplus = 0.5 * dy * sqrt(sum(abs(shear(i,j-1,k,5:8))))/nu_a
                damping(i,j-1,k) = min(damping(i,j-1,k),1.-exp(-(yplus*0.04)**3.))
              endif
              if(lnorm_z(i,j,k)) then
              endif
            endif  !piomelli wall damping
!cstep          endif   !lwallfunc, merged it with k=1 
        end do  !k
      end do    !j
    end do      !i
    !write(6,*) 'Lowest full level reached'
!cstep    if(lwallfunc) then !< special treatment for lowest full level: k=1
      
      do i=2,i1
        do j=2,j1

          if (lnorm_x(i,j,1)) then     !< Wall in x-direction
            emmo = 0.25  * ( &
                ekm(i,j,1)+ekm(i,j-1,1)+ekm(i-1,j-1,1)+ekm(i-1,j,1)  )

            !if(mindist(i,j,1)-0.5*dx<1e-6) write(6,*) 'mindist equal to dx/2'
            call wallaw(v0(i-1,j,1),0.5*dx,nu_a,shear(i-1,j,1,1))
            call wallaw(v0(i,j,1),  0.5*dx,nu_a,shear(i,j,1,2))

            tempvp(i-1,j,1) = tempvp(i-1,j,1) - 0.5 * emmo*((v0(i,j,1)-v0(i-1,j,1))/dx) / dx - 0.5 * shear(i-1,j,1,1)/dx
            tempvp(i,j,1)   = tempvp(i,j,1)   + 0.5 * emmo*((v0(i,j,1)-v0(i-1,j,1))/dx) / dx - 0.5 * shear(i,j,1,2)  /dx
            !vp(i-1,j,1) = vp(i-1,j,1) - 0.5 * emmo*((v0(i,j,1)-v0(i-1,j,1))/dx) / dx - 0.5 * shear(i-1,j,1,1)/dx
            !vp(i,j,1)   = vp(i,j,1)   + 0.5 * emmo*((v0(i,j,1)-v0(i-1,j,1))/dx) / dx - 0.5 * shear(i,j,1,2)  /dx

            empo = 0.25  * ( &
                ekm(i,j+1,1)+ekm(i,j,1)+ekm(i-1,j,1)+ekm(i-1,j+1,1)  )

            !if(mindist(i,j,1)-0.5*dx<1e-6) write(6,*) 'mindist equal to dx/2'
            call wallaw(v0(i-1,j+1,1),0.5*dx,nu_a,shear(i-1,j+1,1,1))
            call wallaw(v0(i,j+1,1),  0.5*dx,nu_a,shear(i,j+1,1,2))

            tempvp(i-1,j+1,1) = tempvp(i-1,j+1,1) - 0.5 * emmo*((v0(i,j+1,1)-v0(i-1,j+1,1))/dx) / dx - 0.5 * shear(i-1,j+1,1,1)/dx
            tempvp(i,j+1,1)   = tempvp(i,j+1,1)   + 0.5 * emmo*((v0(i,j+1,1)-v0(i-1,j+1,1))/dx) / dx - 0.5 * shear(i,j+1,1,2)  /dx
            !vp(i-1,j+1,1) = vp(i-1,j+1,1) - 0.5 * emmo*((v0(i,j+1,1)-v0(i-1,j+1,1))/dx) / dx - 0.5 * shear(i-1,j+1,1,1)/dx
            !vp(i,j+1,1)   = vp(i,j+1,1)   + 0.5 * emmo*((v0(i,j+1,1)-v0(i-1,j+1,1))/dx) / dx - 0.5 * shear(i,j+1,1,2)  /dx


            !call xwallscalar(i,j,1,thl0,tempthlp)
            call xwalltemp(i,j,1,thlm,tempthlp)
            do nc=1,nsv
              call xwallscalar(i,j,1,sv0(:,:,:,nc),tempsvp(:,:,:,nc))
            end do
            call xwalle12(i,j,1)
          endif

          if (lnorm_y(i,j,1)) then     !< Wall in y-direction
            emmo = 0.25  * ( &
              ekm(i,j,1)+ekm(i,j-1,1)+ekm(i-1,j-1,1)+ekm(i-1,j,1)  )

            call wallaw(u0(i,j-1,1),0.5*dy,nu_a,shear(i,j-1,1,5))
            call wallaw(u0(i,j,1)  ,0.5*dy,nu_a,shear(i,j,1,6))

            tempup(i,j-1,1) = tempup(i,j-1,1) - 0.5 * emmo * ((u0(i,j,1)-u0(i,j-1,1))/dy)/dy - 0.5 * shear(i,j-1,1,5)/dy
            tempup(i,j,1)   = tempup(i,j,1)   + 0.5 * emmo * ((u0(i,j,1)-u0(i,j-1,1))/dy)/dy - 0.5 * shear(i,j,1,6)  /dy
            !up(i,j-1,1) = up(i,j-1,1) - 0.5 * emmo * ((u0(i,j,1)-u0(i,j-1,1))/dy)/dy - 0.5 * shear(i,j-1,1,5)/dy
            !up(i,j,1)   = up(i,j,1)   + 0.5 * emmo * ((u0(i,j,1)-u0(i,j-1,1))/dy)/dy - 0.5 * shear(i,j,1,6)  /dy

            epmo = 0.25  * ( &
              ekm(i+1,j,1)+ekm(i+1,j-1,1)+ekm(i,j-1,1)+ekm(i,j,1)  )

            call wallaw(u0(i+1,j-1,1),0.5*dy,nu_a,shear(i+1,j-1,1,5))
            call wallaw(u0(i+1,j,1)  ,0.5*dy,nu_a,shear(i+1,j,1,6))

            tempup(i+1,j-1,1) = tempup(i+1,j-1,1) - 0.5 * emmo * ((u0(i+1,j,1)-u0(i+1,j-1,1))/dy)/dy - 0.5 * shear(i+1,j-1,1,5)/dy
            tempup(i+1,j,1)   = tempup(i+1,j,1)   + 0.5 * emmo * ((u0(i+1,j,1)-u0(i+1,j-1,1))/dy)/dy - 0.5 * shear(i+1,j,1,6)  /dy
            !up(i+1,j-1,1) = up(i+1,j-1,1) - 0.5 * emmo * ((u0(i+1,j,1)-u0(i+1,j-1,1))/dy)/dy - 0.5 * shear(i+1,j-1,1,5)/dy
            !up(i+1,j,1)   = up(i+1,j,1)   + 0.5 * emmo * ((u0(i+1,j,1)-u0(i+1,j-1,1))/dy)/dy - 0.5 * shear(i+1,j,1,6)  /dy


           ! call ywallscalar(i,j,1,thl0,tempthlp)
            call ywalltemp(i,j,1,thlm,tempthlp)
            do nc=1,nsv
              call ywallscalar(i,j,1,sv0(:,:,:,nc),tempsvp(:,:,:,nc))
            end do
            call ywalle12(i,j,1)

          endif
          ekmg(i,j)=-ekm(i,j,1)+2.*nu_a
          ! Wall shear with the ground
          ! write(6,*) 'i=',i
          emom = ( dzh(1) * ( ekm(i,j,1)  + ekm(i-1,j,1)  )  + &
                   dzf(1) * ( ekmg(i,j) + ekmg(i-1,j)     ) ) / &
                 ( 4.   * dzh(1) )

          up(i,j,1)   = up(i,j,1)  - shear(i,j,1,10) /dzf(1) + emom * rhobh(1)/rhobf(1) *((u0(i,j,1)-u0g(i,j))/dzh(1))/dzh(1)

          eomm = ( dzh(1) * ( ekm(i,j,1)  + ekm(i,j-1,1)  )  + &
                   dzf(1) * ( ekmg(i,j)   + ekmg(i,j-1)   ) ) / ( 4.  * dzh(1) )

          vp(i,j,1)   = vp(i,j,1) - shear(i,j,1,12)  /dzf(1)   + eomm * rhobh(1)/rhobf(1) *((v0(i,j,1)-v0g(i,j))/dzh(1))/dzh(1)

          yplus = 0.5*dzh(1)*sqrt(sum(abs(shear(i,j,1,9:12))))/nu_a
          damping(i,j,1)   = min(damping(i,j,1),  1.-exp(-(yplus*0.04)**3.))
        end do
      end do
      !if(myid==1) write(6,*) 'up1 in modibm voor lnormx',up(5,2,2)
      !if(myid==0) write(6,*) 'up0 in modibm voor lnormx',up(8,6,2)
      do i=1,i1
        do j=1,j1
          do k=1,kmax
           !libm ipos=i+myidx*imax
           !libm jpos=j+myidy*jmax
           !libm if (limmersed_boundary(ipos,jpos,k)) then
            if (libm(i,j,k)) then
               up(i,j,k)=-um(i,j,k)*rk3coefi
              vp(i,j,k)=-vm(i,j,k)*rk3coefi
              wp(i,j,k)=-wm(i,j,k)*rk3coefi
              thlp(i,j,k)=(thlwall-thlm(i,j,k))*rk3coefi
              e12p(i,j,k)=(e12min-e12m(i,j,k))*rk3coefi
            else
              up(i,j,k)=up(i,j,k)+tempup(i,j,k)
              vp(i,j,k)=vp(i,j,k)+tempvp(i,j,k)
              wp(i,j,k)=wp(i,j,k)+tempwp(i,j,k)
              thlp(i,j,k)=thlp(i,j,k)+tempthlp(i,j,k)
            endif
            !< Now the normal velocities
            if (lnorm_x(i,j,k)) then
              !write(6,*) 'lnorm_x found, um = ',um(i,j,k)
              !if(myid==0 .and. (i==8 .and. j==6) .and. k==2) write(6,*) 'Before: u0 = ',u0(i,j,k),' um = ',um(i,j,k),' up = ',up(i,j,k)
              up(i,j,k)=-um(i,j,k)*rk3coefi
              !if((i==2 .and. j==4) .and. k==1) write(6,*) 'After: u0 = ',u0(i,j,k),' um = ',um(i,j,k),' up = ',up(i,j,k),' rk3coef = ',rk3coef
            end if
            if (lnorm_y(i,j,k)) then
              !write(6,*) 'lnorm_y found, vm = ',vm(i,j,k)
              vp(i,j,k)=-vm(i,j,k)*rk3coefi
            end if
            if (lnorm_z(i,j,k)) then
              !write(6,*) 'lnorm_z found, wm = ',wm(i,j,k)
              !wp(i,j,k)=-wm(i,j,k)*rk3coefi
            end if
          end do
        end do
      end do
    endif  !lwallfunc

    call excjs( up  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( vp  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( wp  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( damping  , 2,i1,2,j1,1,kmax,ih,jh)
    call excjs( e12p  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( thlp  , 2,i1,2,j1,1,k1,ih,jh)
    do nc=1,nsv
      svp(:,:,:,nc)=svp(:,:,:,nc)+tempsvp(:,:,:,nc)
      call excjs( svp(:,:,:,nc)  , 2,i1,2,j1,1,k1,ih,jh)
    enddo


    if(maxval(sqrt(u0**2.+v0**2.+w0**2.))>16) then
      maxlocx=maxloc(u0**2.+v0**2.+w0**2.)
      !write(6,*) 'maxlocx = ',maxlocx(1),maxlocx(2),maxlocx(3)
      !write(6,*) 'ERROR: vel>16, maxloc = ',maxloc(sqrt(u0**2.+v0**2.+w0**2.)), 'u0 and um and up are here:',u0(maxlocx(1),maxlocx(2),maxlocx(3)),um(maxlocx(1),maxlocx(2),maxlocx(3)),up(maxlocx(1),maxlocx(2),maxlocx(3))
      !write(6,*) 'v0 and vm and vp are here:',v0(maxlocx(1),maxlocx(2),maxlocx(3)),vm(maxlocx(1),maxlocx(2),maxlocx(3)),vp(maxlocx(1),maxlocx(2),maxlocx(3))
      !write(6,*) 'w0 and wm and wp are here:',w0(maxlocx(1),maxlocx(2),maxlocx(3)),wm(maxlocx(1),maxlocx(2),maxlocx(3)),wp(maxlocx(1),maxlocx(2),maxlocx(3))
      !write(6,*) 'timee = ',timee
    endif
    return
  end subroutine applyibm


  subroutine zerowallvelocity(simid) !<- MK: Set velocity at the immersed boundaries to 0 for a better interaction with the poissonsolver

    use modfields,      only : um, vm, wm, up, vp, wp
    use modglobal,      only : rk3step, imax, jmax, kmax, i1, j1, k1, ih, jh, rdt
    use modmpi,         only : excjs   !cstep ,myidx,myidy
    use modnudgeboundary, only : Nsim

    implicit none
    integer  :: i, j, k,ipos,jpos
    real     :: rk3coef,rk3coefi
integer, intent(in) :: simid

    if (.not. (lapply_ibm .and. simid==Nsim )) return

    rk3coef = rdt / (4. - dble(rk3step))
    rk3coefi = 1. / rk3coef

    do i=1,i1
      do j=1,j1
        do k=1,kmax
          !libm ipos=i+myidx*imax
          !!libm jpos=j+myidy*jmax
          !libm if (limmersed_boundary(ipos,jpos,k)) then
          if (libm(i,j,k)) then
            up(i,j,k)=-um(i,j,k)*rk3coefi
            vp(i,j,k)=-vm(i,j,k)*rk3coefi
            wp(i,j,k)=-wm(i,j,k)*rk3coefi
          endif
          if (lnorm_x(i,j,k)) then
            up(i,j,k)=-um(i,j,k)*rk3coefi
          end if
          if (lnorm_y(i,j,k)) then
            vp(i,j,k)=-vm(i,j,k)*rk3coefi
          end if
          if (lnorm_z(i,j,k)) then
          ! Not needed here
          end if
        end do
      end do
    end do
    call excjs( up  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( vp  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( wp  , 2,i1,2,j1,1,k1,ih,jh)

    return
  end subroutine zerowallvelocity


  subroutine xwallscalar(i,j,k,putin,putout)

    use modglobal,      only : ih, i1, jh, j1, k1, dx2i
    use modsubgriddata, only : ekh
    use modfields, only : sv0
!cstep    use modmpi, only : myid

    implicit none

    integer, intent(in)    :: i,j,k
    real,    intent(in)    :: putin(2-ih:i1+ih,2-jh:j1+jh,k1)
    real,    intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)

    ! zero flux through the boundary is applied here
    if(.not. (putin(i,j,k)==0)) then
      putout(i,j,k)   = putout(i,j,k) &
                      + 0.5 * (ekh(i,j,k)+ekh(i-1,j,k))*(putin(i,j,k)-putin(i-1,j,k))*dx2i
    elseif(.not. (putin(i-1,j,k)==0)) then
      putout(i-1,j,k) = putout(i-1,j,k) &
                      - 0.5 * (ekh(i,j,k)+ekh(i-1,j,k))*(putin(i,j,k)-putin(i-1,j,k))*dx2i
    endif

    return
  end subroutine xwallscalar

  subroutine ywallscalar(i,j,k,putin,putout)

    use modglobal,      only : ih, i1, jh, j1, k1, dy2i
    use modsubgriddata, only : ekh

    implicit none

    integer, intent(in)    :: i,j,k
    real,    intent(in)    :: putin(2-ih:i1+ih,2-jh:j1+jh,k1)
    real,    intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)

    ! zero flux through the boundary is applied here
    if(.not. (putin(i,j,k)==0)) then
      putout(i,j,k)   = putout(i,j,k) &
                      + 0.5 * (ekh(i,j,k)+ekh(i,j-1,k)) *(putin(i,j,k)-putin(i,j-1,k))*dy2i
    elseif(.not. (putin(i,j-1,k)==0)) then
      putout(i,j-1,k) = putout(i,j-1,k) &
                      - 0.5 * (ekh(i,j,k)+ekh(i,j-1,k)) *(putin(i,j,k)-putin(i,j-1,k))*dy2i
    endif

    return
  end subroutine ywallscalar

  subroutine zwallscalar(i,j,k,putin,putout)

    use modglobal,      only : ih, i1, jh, j1, k1, dzf, dzh
    use modsubgriddata, only : ekh
    use modfields,      only : rhobh, rhobf, sv0

    implicit none

    integer, intent(in)    :: i,j,k
    real,    intent(in)    :: putin(2-ih:i1+ih,2-jh:j1+jh,k1)
    real,    intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)

    ! zero flux through the boundary is applied here
    if(.not. (putin(i,j,k)==0)) then
      putout(i,j,k)   = putout(i,j,k) &
                      + 0.5 * rhobh(k)/rhobf(k-1) * (dzf(k-1)*ekh(i,j,k) + dzf(k)*ekh(i,j,k-1)) &
                      *  (putin(i,j,k)-putin(i,j,k-1)) / dzh(k)**2 /dzf(k-1)
    elseif(.not. (putin(i,j,k-1)==0)) then
      putout(i,j,k-1) = putout(i,j,k-1) &
                      - 0.5 * rhobh(k)/rhobf(k-1) * (dzf(k-1)*ekh(i,j,k) + dzf(k)*ekh(i,j,k-1)) &
                      *  (putin(i,j,k)-putin(i,j,k-1)) / dzh(k)**2 /dzf(k-1)
    endif

    return
  end subroutine zwallscalar

  subroutine xwalltemp(i,j,k,thlm,tempthlp)

    use modglobal,      only : ih, i1, jh, j1, k1, rk3step, rdt

    implicit none

    integer, intent(in)    :: i,j,k
    real,    intent(in)    :: thlm(2-ih:i1+ih,2-jh:j1+jh,k1)
    real,    intent(inout) :: tempthlp(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                   :: rk3coef,rk3coefi

    rk3coef = rdt / (4. - dble(rk3step))
    rk3coefi = 1. / rk3coef

    tempthlp(i  ,j,k) = ct * (thlwall - thlm(i  ,j,k)) * rk3coefi
    tempthlp(i-1,j,k) = ct * (thlwall - thlm(i-1,j,k)) * rk3coefi

    return
  end subroutine xwalltemp

  subroutine ywalltemp(i,j,k,thlm,tempthlp)

    use modglobal,      only : ih, i1, jh, j1, k1, rk3step, rdt

    implicit none

    integer, intent(in)    :: i,j,k
    real,    intent(in)    :: thlm(2-ih:i1+ih,2-jh:j1+jh,k1)
    real,    intent(inout) :: tempthlp(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                   :: rk3coef,rk3coefi

    rk3coef = rdt / (4. - dble(rk3step))
    rk3coefi = 1. / rk3coef

    tempthlp(i,j  ,k) = ct * (thlwall - thlm(i,j  ,k)) * rk3coefi
    tempthlp(i,j-1,k) = ct * (thlwall - thlm(i,j-1,k)) * rk3coefi

    return
  end subroutine ywalltemp


  subroutine xwalle12(i,j,k)

    use modglobal,      only : ih, i1, jh, j1, k1, dx2i, dx, dy, dzh
    use modsubgriddata, only : ekm
    use modfields,      only : e12p, e120, u0, v0, w0

    implicit none

    integer, intent(in)    :: i,j,k

    if(.not.(k==1)) then
      if(.not. (e12p(i,j,k)==0)) then
        e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
                                  + ekm(i,j,k)/(2*e120(i,j,k))* (&  !source terms
                                     -((w0(i,j,k+1)-w0(i-1,j,k+1))  / dx             + &
                                       (u0(i,j,k+1)-u0(i,j,k))      / dzh(k+1) )**2  + &
                                     +(2.*(w0(i,j,k+1))             / dx             + &
                                       (u0(i,j,k+1)-u0(i,j,k))      / dzh(k+1) )**2  + &

                                     -((w0(i,j,k)-w0(i-1,j,k))      / dx             + &
                                       (u0(i,j,k)-u0(i,j,k-1))      / dzh(k)   )**2  + &
                                     +(2.*(w0(i,j,k))               / dx             + &
                                       (u0(i,j,k)-u0(i,j,k-1))      / dzh(k)   )**2  + &

                                     -((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
                                       (v0(i,j+1,k)-v0(i-1,j+1,k))  / dx       )**2  + &
                                     +((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
                                       (2.*v0(i,j+1,k))             / dx       )**2  + &

                                     -((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
                                       (v0(i,j,k)-v0(i-1,j,k))      / dx       )**2  + &
                                     +((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
                                       (2.*v0(i,j,k))               / dx       )**2    &
                                    )
      elseif(.not. (e12p(i-1,j,k)==0)) then
        e12p(i-1,j,k) = e12p(i-1,j,k) + (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
                                  + ekm(i-1,j,k)/(2*e120(i-1,j,k))* (&  !source terms
                                       -((w0(i,j,k)-w0(i-1,j,k))    / dx             + &
                                         (u0(i,j,k)-u0(i,j,k-1))    / dzh(k)   )**2  + &
                                       +((-2.*w0(i-1,j,k))          / dx             + &
                                         (u0(i,j,k)-u0(i,j,k-1))    / dzh(k)   )**2  + &

                                       -((w0(i,j,k+1)-w0(i-1,j,k+1))/ dx             + &
                                         (u0(i,j,k+1)-u0(i,j,k))    / dzh(k+1) )**2  + &
                                       +((-2.*w0(i-1,j,k+1))        / dx             + &
                                         (u0(i,j,k+1)-u0(i,j,k))    / dzh(k+1) )**2  + &

                                       -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
                                       +((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (-2.*v0(i-1,j,k))          / dx       )**2  + &

                                       -((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
                                         (v0(i,j+1,k)-v0(i-1,j+1,k))/ dx       )**2  + &
                                       +((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
                                         (-2.*v0(i-1,j+1,k))        / dx       )**2    &
                                    )
    endif
    else !Special treatment for the lowest full level: k=1
      if(.not. (e12p(i,j,k)==0)) then
        e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
                                  + ekm(i,j,k)/(2*e120(i,j,k))* (&  !source terms
                                     -((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
                                       (v0(i,j+1,k)-v0(i-1,j+1,k))  / dx       )**2  + &
                                     +((u0(i,j+1,k)-u0(i,j,k))      / dy             + &
                                       (2.*v0(i,j+1,k))             / dx       )**2  + &

                                     -((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
                                       (v0(i,j,k)-v0(i-1,j,k))      / dx       )**2  + &
                                     +((u0(i,j,k)-u0(i,j-1,k))      / dy             + &
                                       (2.*v0(i,j,k))               / dx       )**2    &
                                    )
      elseif(.not. (e12p(i-1,j,k)==0)) then
        e12p(i-1,j,k) = e12p(i-1,j,k) + (ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k))*dx2i &
                                  + ekm(i-1,j,k)/(2*e120(i-1,j,k))* (&  !source terms

                                       -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
                                       +((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (-2.*v0(i-1,j,k))          / dx       )**2  + &

                                       -((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
                                         (v0(i,j+1,k)-v0(i-1,j+1,k))/ dx       )**2  + &
                                       +((u0(i,j+1,k)-u0(i,j,k))    / dy             + &
                                         (-2.*v0(i-1,j+1,k))        / dx       )**2    &
                                    )
      endif
    endif
  end subroutine xwalle12


  subroutine ywalle12(i,j,k)

    use modglobal,      only : ih, i1, jh, j1, k1, dy2i, dx, dy, dzh
    use modsubgriddata, only : ekm
    use modfields,      only : e12p, e120, u0, v0, w0

    implicit none

    integer, intent(in)    :: i,j,k

    if(.not.(k==1)) then
      if(.not. (e12p(i,j,k)==0)) then
        e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
                                  + ekm(i,j,k)/(2.*e120(i,j,k))* (&  !source terms
                                       -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
                                       +((2.*u0(i,j,k))             / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

                                       -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
                                       +((2.*u0(i+1,j,k))           / dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &

                                       -((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
                                         (w0(i,j,k+1)-w0(i,j-1,k+1))/ dy       )**2  + &
                                       +((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
                                         (2.*w0(i,j,k+1))           / dy       )**2  + &

                                       -((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
                                         (w0(i,j,k)-w0(i,j-1,k))    / dy       )**2  + &
                                       +((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
                                         (2.*w0(i,j,k))             / dy       )**2    &
                                    )
      elseif(.not. (e12p(i,j-1,k)==0)) then
        e12p(i,j-1,k) = e12p(i,j-1,k) + (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
                                  + ekm(i,j-1,k)/(2.*e120(i,j-1,k))* (&  !source terms
                                       -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
                                       +((-2.*u0(i,j-1,k))          / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

                                       -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
                                       +((-2.*u0(i+1,j-1,k))        / dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &

                                       -((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
                                         (w0(i,j,k)-w0(i,j-1,k))    / dy       )**2  + &
                                       +((v0(i,j,k)-v0(i,j,k-1))    / dzh(k)         + &
                                         (-2.*w0(i,j-1,k))          / dy       )**2  + &

                                       -((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
                                         (w0(i,j,k+1)-w0(i,j-1,k+1))/ dy       )**2  + &
                                       +((v0(i,j,k+1)-v0(i,j,k))    / dzh(k+1)       + &
                                         (-2.*w0(i,j-1,k+1))        / dy       )**2    &
                                    )
      endif
    else !Special treatment for the lowest full level: k=1
      if(.not. (e12p(i,j,k)==0)) then
        e12p(i,j,k)   = e12p(i,j,k)   - (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
                                  + ekm(i,j,k)/(2.*e120(i,j,k))* (&  !source terms
                                       -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
                                       +((2.*u0(i,j,k))             / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

                                       -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
                                       +((2.*u0(i+1,j,k))           / dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2    &
                                    )

      elseif(.not. (e12p(i,j-1,k)==0)) then
        e12p(i,j-1,k) = e12p(i,j-1,k) + (ekm(i,j,k)+ekm(i,j-1,k))*(e120(i,j,k)-e120(i,j-1,k))*dy2i &
                                  + ekm(i,j-1,k)/(2.*e120(i,j-1,k))* (&  !source terms
                                       -((u0(i,j,k)-u0(i,j-1,k))    / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &
                                       +((-2.*u0(i,j-1,k))          / dy             + &
                                         (v0(i,j,k)-v0(i-1,j,k))    / dx       )**2  + &

                                       -((u0(i+1,j,k)-u0(i+1,j-1,k))/ dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2  + &
                                       +((-2.*u0(i+1,j-1,k))        / dy             + &
                                         (v0(i+1,j,k)-v0(i,j,k))    / dx       )**2    &
                                    )
      endif
    endif
  end subroutine ywalle12

   subroutine zwalle12(i,j,k)

    use modglobal,      only : ih, i1, jh, j1, k1, dzf, dx, dy, dzh
    use modsubgriddata, only : ekm
    use modfields,      only : e12p, e120, u0, v0, w0, rhobh, rhobf

    implicit none

    integer, intent(in)    :: i,j,k

    if(.not. (e12p(i,j,k)==0)) then
    e12p(i,j,k)   = e12p(i,j,k) - rhobh(k)/rhobf(k-1) * (dzf(k)*ekm(i,j,k-1) + dzf(k-1)*ekm(i,j,k)) &
                                  *(e120(i,j,k)-e120(i,j,k-1)) / dzh(k)**2 /dzf(k)     &
                                + ekm(i,j,k)/(2.*e120(i,j,k))* (&  !source terms
                                   -((w0(i,j,k)-w0(i-1,j,k))        / dx             + &
                                     (u0(i,j,k)-u0(i,j,k-1))        / dzh(k)   )**2  + &
                                   +((w0(i,j,k)-w0(i-1,j,k))        / dx             + &
                                     (2.*u0(i,j,k))                 / dzh(k)   )**2  + &

                                   -((w0(i+1,j,k)-w0(i,j,k))        / dx             + &
                                     (u0(i+1,j,k)-u0(i+1,j,k-1))    / dzh(k)   )**2  + &
                                   +((w0(i+1,j,k)-w0(i,j,k))        / dx             + &
                                     (2.*u0(i+1,j,k))               / dzh(k)   )**2  + &

                                   -((v0(i,j,k)-v0(i,j,k-1))        / dzh(k)         + &
                                     (w0(i,j,k)-w0(i,j-1,k))        / dy       )**2  + &
                                   +((2.*v0(i,j,k))                 / dzh(k)         + &
                                     (w0(i,j,k)-w0(i,j-1,k))        / dy       )**2  + &

                                   -((v0(i,j+1,k)-v0(i,j+1,k-1))    / dzh(k)         + &
                                     (w0(i,j+1,k)-w0(i,j,k))        / dy       )**2  + &
                                   +(2.*v0(i,j+1,k)                 / dzh(k)         + &
                                     (w0(i,j+1,k)-w0(i,j,k))        / dy       )**2    &
                                )
    elseif(.not. (e12p(i,j,k-1)==0)) then
      e12p(i,j,k-1) = e12p(i,j,k-1) + rhobh(k)/rhobf(k) * (dzf(k-1)*ekm(i,j,k) + dzf(k)*ekm(i,j,k-1)) &
                                    *(e120(i,j,k)-e120(i,j,k-1)) / dzh(k)**2  /dzf(k)  &
                                  + ekm(i,j,k-1)/(2.*e120(i,j,k-1))* (&  !source terms
                                   -((w0(i,j,k)-w0(i-1,j,k))        / dx             + &
                                     (u0(i,j,k)-u0(i,j,k-1))        / dzh(k)   )**2  + &
                                   +((w0(i,j,k)-w0(i-1,j,k))        / dx             + &
                                     (-2.*u0(i,j,k-1))              / dzh(k)   )**2  + &

                                   -((w0(i+1,j,k)-w0(i,j,k))        / dx             + &
                                     (u0(i+1,j,k)-u0(i+1,j,k-1))    / dzh(k)   )**2  + &
                                   +((w0(i+1,j,k)-w0(i,j,k))        / dx             + &
                                     (-2.*u0(i+1,j,k-1))            / dzh(k)   )**2  + &
                                   -((v0(i,j,k)-v0(i,j,k-1))        / dzh(k)         + &
                                     (w0(i,j,k)-w0(i,j-1,k))        / dy       )**2  + &
                                   +((-2.*v0(i,j,k-1))              / dzh(k)         + &
                                     (w0(i,j,k)-w0(i,j-1,k))        / dy       )**2  + &

                                   -((v0(i,j+1,k)-v0(i,j+1,k-1))    / dzh(k)         + &
                                     (w0(i,j+1,k)-w0(i,j,k))        / dy       )**2  + &
                                   +((-2.*v0(i,j+1,k-1))            / dzh(k)         + &
                                     (w0(i,j+1,k)-w0(i,j,k))        / dy       )**2    &
                                )
    endif
  end subroutine zwalle12

  subroutine wallaw(utan,dx,visc,tau)  !< Copied from work by Jasper Tomas
    !use modglobal,         only :

    implicit none

    real, intent(in)  :: utan  ! tangential velocity component
    real, intent(in)  :: dx    ! distance to the wall
    real, intent(in)  :: visc  ! viscosity near the wall
    real, intent(out) :: tau   !

    real dxi,dx5
    real utanabs, utankr, dutan, sub
    real tausub, taupow

    real const1, const2, const3, const4
    real aaa, bbb

    parameter(aaa = 8.3)
    parameter(bbb = 0.1428571429)

    !write(6,*) 'Starting with wallaw'

    dxi = 1./dx
    dx5 = 0.5*dx

    const1 = 0.5 * (1. - bbb) * aaa ** ((1. + bbb) / (1. - bbb))
    const2 = (1. + bbb) / aaa
    const3 = aaa ** (2. / (1. - bbb))
    const4 = 2. / (1. + bbb)

    utanabs = abs(utan)
    utankr  = 0.5 * visc * dxi * const3
    dutan   = utankr - utanabs
    sub     = max (sign(1.,dutan),0.) ! sub = 1 for viscous sublayer and 0 for outside.

    tausub  = 2. * visc * utanabs * dxi
    taupow  = (const1 * (visc * dxi) ** (1. + bbb) + const2 * ((visc * dxi) ** bbb) * utanabs) ** const4

    tau = sub * tausub + (1. - sub) * taupow ! use linear profile inside viscous sublayer and use power law outside
    tau = sign(tau,utan)   ! give tau the same sign as utan
    return
  end subroutine wallaw

end module modibm
