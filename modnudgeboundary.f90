!> \file modnudgeboundary.f90
!! By Pim van Dorp (PVD), TU Delft, section Atmospheric Physics, 10 dec 2015
!! Nudge boundary to prescribed values

module modnudgeboundary 
  use modglobal, only : longint

  implicit none 
  logical :: lnudgeboundary = .false. 
  logical :: lstatref = .false. 

  integer :: Nsim = 1 
  integer :: statid = 1
  integer :: refid = 1
  integer :: turid = 1     

  real, allocatable :: fnudgeglob(:,:,:) ! global array of fnudge values
  real, allocatable :: fnudgeloc(:,:,:) ! local, cpu dependent array of fnudge values

  integer :: nudgedepthgr = 10 ! number of nudge grid points
  
  ! Prognostic variables; first dimension: 1=turbine simulation, 2=reference simulation
  real, allocatable :: umsave(:,:,:,:)        !<   x-component of velocity at time step t-1
  real, allocatable :: vmsave(:,:,:,:)        !<   y-component of velocity at time step t-1
  real, allocatable :: wmsave(:,:,:,:)        !<   z-component of velocity at time step t-1
  real, allocatable :: thlmsave(:,:,:,:)      !<   liq. water pot. temperature at time step t-1
  real, allocatable :: e12msave(:,:,:,:)      !<   square root of turb. kin. energy at time step t-1
  real, allocatable :: qtmsave(:,:,:,:)       !<   total specific humidity at time step t

  real, allocatable :: u0save(:,:,:,:)        !<   x-component of velocity at time step t
  real, allocatable :: v0save(:,:,:,:)        !<   y-component of velocity at time step t
  real, allocatable :: w0save(:,:,:,:)        !<   z-component of velocity at time step t
  real, allocatable :: thl0save(:,:,:,:)      !<   liq. water pot. temperature at time step t
  real, allocatable :: qt0save(:,:,:,:)       !<   total specific humidity at time step t
  real, allocatable :: ql0save(:,:,:,:)   
  real, allocatable :: ql0hsave(:,:,:,:)  
  real, allocatable :: e120save(:,:,:,:)      !<   square root of turb. kin. energy at time step t

  real, allocatable :: dthvdzsave(:,:,:,:)  
  real, allocatable :: ekmsave(:,:,:,:)  
  real, allocatable :: tmp0save(:,:,:,:)  
  real, allocatable :: eslsave(:,:,:,:)  
  real, allocatable :: qvslsave(:,:,:,:)  
  real, allocatable :: qvsisave(:,:,:,:)  

  real, allocatable :: thv0hsave(:,:,:,:)  

  real, allocatable :: presfsave(:,:)  
  real, allocatable :: preshsave(:,:)  

  real, allocatable :: thvhsave(:,:)  

  real, allocatable :: u0avsave(:,:)
  real, allocatable :: v0avsave(:,:)
  real, allocatable :: thl0avsave(:,:)
  real, allocatable :: qt0avsave(:,:)

  integer, allocatable :: ksfcsave(:,:,:)    !< !cstep ground surface level, needed for obstacles with immersed boundary method
!cstep implements this or the other ibm booleans?  logical, allocatable :: libmsave(:,:,:)    

contains
  subroutine initnudgeboundary
    use modglobal, only: itot, jtot, kmax, i1, i2, j1, j2, k1, ih, jh,ifnamopt, fname_options, tres, ladaptive, dtmax,btime, dx, dy , pi, dt
    use modmpi, only : myidx,myidy,myid,MPI_INTEGER,MY_REAL, MPI_SUM,MPI_COMM_WORLD,MPI_LOGICAL,comm3d,mpierr
    use modfields, only : u0, v0, w0, e120, thl0, qt0, ql0, ql0h, tmp0, &
                          um, vm, wm, e12m, thlm, qtm, &
                          presf, presh, dthvdz, esl, qvsl, qvsi, thv0h, &
                          u0av, v0av, thl0av, qt0av, thvh, & !cstep IBM
                          ksfc
    use modsubgrid, only : ekm

    implicit none
    integer i,j,k,ierr, simid


    namelist/NUDGEBOUNDARY/ lnudgeboundary, lstatref, nudgedepthgr

    if (myid==0) then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
        read (ifnamopt,NUDGEBOUNDARY,iostat=ierr)
        if (ierr > 0) then
          print *, 'Problem in namoptions NUDGEBOUNDARY'
          print *, 'iostat error: ', ierr
          stop 'ERROR: Problem in namoptions NUDGEBOUNDARY'
        endif
        write(6 ,NUDGEBOUNDARY)
      close(ifnamopt)
    end if

    call MPI_BCAST(lnudgeboundary ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr) 
    call MPI_BCAST(lstatref,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr) 
    call MPI_BCAST(nudgedepthgr       ,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr) 

    if (.not. (lnudgeboundary)) return

    if (lnudgeboundary) Nsim = 2
    if (lnudgeboundary) turid = 2
    if (.not. lstatref .and. Nsim == 2) statid = 2

    allocate(umsave(2-ih:i1+ih,2-jh:j1+jh,k1,2))
    allocate(vmsave(2-ih:i1+ih,2-jh:j1+jh,k1,2))
    allocate(wmsave(2-ih:i1+ih,2-jh:j1+jh,k1,2))
    allocate(thlmsave(2-ih:i1+ih,2-jh:j1+jh,k1,2))
    allocate(e12msave(2-ih:i1+ih,2-jh:j1+jh,k1,2))
    allocate(qtmsave(2-ih:i1+ih,2-jh:j1+jh,k1,2))

    allocate(u0save(2-ih:i1+ih,2-jh:j1+jh,k1,2))
    allocate(v0save(2-ih:i1+ih,2-jh:j1+jh,k1,2))
    allocate(w0save(2-ih:i1+ih,2-jh:j1+jh,k1,2))
    allocate(thl0save(2-ih:i1+ih,2-jh:j1+jh,k1,2))
    allocate(qt0save(2-ih:i1+ih,2-jh:j1+jh,k1,2))
    allocate(ql0save(2-ih:i1+ih,2-jh:j1+jh,k1,2))
    allocate(ql0hsave(2-ih:i1+ih,2-jh:j1+jh,k1,2))
    allocate(e120save(2-ih:i1+ih,2-jh:j1+jh,k1,2))
    allocate(dthvdzsave(2-ih:i1+ih,2-jh:j1+jh,k1,2))
    allocate(ekmsave(2-ih:i1+ih,2-jh:j1+jh,k1,2))
    allocate(tmp0save(2-ih:i1+ih,2-jh:j1+jh,k1,2))
    allocate(eslsave(2-ih:i1+ih,2-jh:j1+jh,k1,2))
    allocate(qvslsave(2-ih:i1+ih,2-jh:j1+jh,k1,2))
    allocate(qvsisave(2-ih:i1+ih,2-jh:j1+jh,k1,2))

    allocate(thv0hsave(2-ih:i1+ih,2-jh:j1+jh,k1,2))

    allocate(presfsave(k1,2))
    allocate(preshsave(k1,2))

    allocate(thvhsave(k1,2))

    allocate(u0avsave(k1,2))
    allocate(v0avsave(k1,2))
    allocate(thl0avsave(k1,2))
    allocate(qt0avsave(k1,2))
   
    allocate(ksfcsave(2-ih:i1+ih,2-jh:j1+jh,2)) 

    do simid=1,2
      umsave(:,:,:,simid) = um(:,:,:)
      vmsave(:,:,:,simid) = vm(:,:,:)
      wmsave(:,:,:,simid) = wm(:,:,:)
      e12msave(:,:,:,simid) = e12m(:,:,:)
      thlmsave(:,:,:,simid) = thlm(:,:,:)
      qtmsave(:,:,:,simid) = qtm(:,:,:)

      u0save(:,:,:,simid) = u0(:,:,:)
      v0save(:,:,:,simid) = v0(:,:,:)
      w0save(:,:,:,simid) = w0(:,:,:)
      thl0save(:,:,:,simid) = thl0(:,:,:)
      qt0save(:,:,:,simid) = qt0(:,:,:)
      ql0save(:,:,:,simid) = ql0(:,:,:)
      ql0hsave(:,:,:,simid) = ql0h(:,:,:)
      e120save(:,:,:,simid) = e120(:,:,:)
      dthvdzsave(:,:,:,simid) = dthvdz(:,:,:)
      ekmsave(:,:,:,simid) = ekm(:,:,:)
      tmp0save(:,:,:,simid) = tmp0(:,:,:)
      eslsave(:,:,:,simid) = esl(:,:,:)
      qvslsave(:,:,:,simid) = qvsl(:,:,:)
      qvsisave(:,:,:,simid) = qvsi(:,:,:)

      thv0hsave(:,:,:,simid) = thv0h(:,:,:)

      presfsave(:,simid) = presf(:)
      preshsave(:,simid) = presf(:)

      thvhsave(:,simid) = thvh(:)

      u0avsave(:,simid) = u0av(:)
      v0avsave(:,simid) = v0av(:)
      thl0avsave(:,simid) = thl0av(:)
      qt0avsave(:,simid) = qt0av(:)
    end do
    ksfcsave(:,:,1) = 1
    ksfcsave(:,:,2) = ksfc(:,:)

    allocate(fnudgeglob(1-ih:itot+ih,1-jh:jtot+jh,1:k1))
    allocate(fnudgeloc(2-ih:i1+ih,2-jh:j1+jh,k1))

    call calcfnudge

    if (myid == 0 ) then
      write(*,*) 'nudgedepthgr = ', nudgedepthgr 
      write(*,*) 'Succesfully initialized modnudgeboundary'
    end if

  end subroutine initnudgeboundary

  subroutine calcfnudge
    use modglobal, only : pi, itot, jtot, ih, jh, k1, j1, i1, kmax
    use modfields, only : u0av, v0av
!cstep    use modwindturbinedata, only : turhzgr
    use modmpi, only : myidx, myidy

    implicit none
    integer i,j,k
    real fnudge

    fnudgeglob = 0.
    fnudgeloc = 0.

    do i=1,nudgedepthgr 

      fnudge = 0.5 + 0.5*COS((pi/(nudgedepthgr-1))*(i-1))

      fnudgeglob(i,i:jtot-i+1,:) = fnudge
      fnudgeglob(itot-i+1,i:jtot-i+1,:) = fnudge
      fnudgeglob(i+1:(itot-i),i,:) = fnudge
      fnudgeglob(i+1:(itot-i),jtot-i+1,:) = fnudge

    end do

    do k=1,kmax
      do j=2,j1
        do i=2,i1
          fnudgeloc(i,j,k) = fnudgeglob(iglob(i,myidx),jglob(j,myidy),k)
        end do
      end do
    end do

  end subroutine calcfnudge

  subroutine loadfields(simid)
    use modfields, only : u0, v0, w0, e120, thl0, qt0, ql0, ql0h, tmp0, &
                          um, vm, wm, e12m, thlm, qtm, &
                          presf, presh, dthvdz, esl, qvsl, qvsi, thv0h, &
                          u0av, v0av, thl0av, qt0av, thvh, &  !cstep IBM
                          ksfc
    use modsubgrid, only : ekm

    implicit none
    integer, intent(in) :: simid  

    um(:,:,:) = umsave(:,:,:,simid)
    vm(:,:,:) = vmsave(:,:,:,simid)
    wm(:,:,:) = wmsave(:,:,:,simid)
    e12m(:,:,:) = e12msave(:,:,:,simid)
    thlm(:,:,:) = thlmsave(:,:,:,simid)
    qtm(:,:,:) = qtmsave(:,:,:,simid)

    u0(:,:,:) = u0save(:,:,:,simid)
    v0(:,:,:) = v0save(:,:,:,simid)
    w0(:,:,:) = w0save(:,:,:,simid)
    thl0(:,:,:) = thl0save(:,:,:,simid)
    qt0(:,:,:) = qt0save(:,:,:,simid)
    ql0(:,:,:) = ql0save(:,:,:,simid)
    ql0h(:,:,:) = ql0hsave(:,:,:,simid)
    e120(:,:,:) = e120save(:,:,:,simid)
    dthvdz(:,:,:) = dthvdzsave(:,:,:,simid)
    ekm(:,:,:) = ekmsave(:,:,:,simid)
    tmp0(:,:,:) = tmp0save(:,:,:,simid)
    esl(:,:,:) = eslsave(:,:,:,simid)
    qvsl(:,:,:) = qvslsave(:,:,:,simid)
    qvsi(:,:,:) = qvsisave(:,:,:,simid)

    thv0h(:,:,:) = thv0hsave(:,:,:,simid)

    presf(:) = presfsave(:,simid)
    presh(:) = presfsave(:,simid)

    thvh(:) = thvhsave(:,simid)

    u0av(:) = u0avsave(:,simid)
    v0av(:) = v0avsave(:,simid)
    thl0av(:) = thl0avsave(:,simid)
    qt0av(:) = qt0avsave(:,simid)

    ksfc(:,:) = ksfcsave(:,:,simid)

  end subroutine loadfields

  subroutine savefields(simid)
    use modfields, only : u0, v0, w0, e120, thl0, qt0, ql0, ql0h, tmp0, &
                          um, vm, wm, e12m, thlm, qtm, &
                          presf, presh, dthvdz, esl, qvsl, qvsi, thv0h, &
                          u0av, v0av, thl0av, qt0av, thvh
    use modsubgrid, only : ekm

    implicit none
    integer, intent(in) :: simid  

    umsave(:,:,:,simid) = um(:,:,:)
    vmsave(:,:,:,simid) = vm(:,:,:)
    wmsave(:,:,:,simid) = wm(:,:,:)
    e12msave(:,:,:,simid) = e12m(:,:,:)
    thlmsave(:,:,:,simid) = thlm(:,:,:)
    qtmsave(:,:,:,simid) = qtm(:,:,:)

    u0save(:,:,:,simid) = u0(:,:,:)
    v0save(:,:,:,simid) = v0(:,:,:)
    w0save(:,:,:,simid) = w0(:,:,:)
    thl0save(:,:,:,simid) = thl0(:,:,:)
    qt0save(:,:,:,simid) = qt0(:,:,:)
    ql0save(:,:,:,simid) = ql0(:,:,:)
    ql0hsave(:,:,:,simid) = ql0h(:,:,:)
    e120save(:,:,:,simid) = e120(:,:,:)
    dthvdzsave(:,:,:,simid) = dthvdz(:,:,:)
    ekmsave(:,:,:,simid) = ekm(:,:,:)
    tmp0save(:,:,:,simid) = tmp0(:,:,:)
    eslsave(:,:,:,simid) = esl(:,:,:)
    qvslsave(:,:,:,simid) = qvsl(:,:,:)
    qvsisave(:,:,:,simid) = qvsi(:,:,:)

    thv0hsave(:,:,:,simid) = thv0h(:,:,:)

    presfsave(:,simid) = presf(:)
    preshsave(:,simid) = presf(:)

    thvhsave(:,simid) = thvh(:)

    u0avsave(:,simid) = u0av(:)
    v0avsave(:,simid) = v0av(:)
    thl0avsave(:,simid) = thl0av(:)
    qt0avsave(:,simid) = qt0av(:)

  end subroutine savefields

  subroutine nudgeboundary 
    use modglobal, only : kmax, i1, j1, rdt
    use modfields, only : u0, v0, w0, thl0, e120, qt0, up, vp, wp, thlp, qtp, e12p

    implicit none
    integer i,j,k

    do k=1,kmax
      do j=2,j1
        do i=2,i1
          up(i,j,k) = (1-fnudgeloc(i,j,k))*up(i,j,k) + fnudgeloc(i,j,k)*(u0save(i,j,k,1)-u0(i,j,k))/rdt
          vp(i,j,k) = (1-fnudgeloc(i,j,k))*vp(i,j,k) + fnudgeloc(i,j,k)*(v0save(i,j,k,1)-v0(i,j,k))/rdt
          wp(i,j,k) = (1-fnudgeloc(i,j,k))*wp(i,j,k) + fnudgeloc(i,j,k)*(w0save(i,j,k,1)-w0(i,j,k))/rdt
          thlp(i,j,k) = (1-fnudgeloc(i,j,k))*thlp(i,j,k) + fnudgeloc(i,j,k)*(thl0save(i,j,k,1)-thl0(i,j,k))/rdt
          qtp(i,j,k) = (1-fnudgeloc(i,j,k))*qtp(i,j,k) + fnudgeloc(i,j,k)*(qt0save(i,j,k,1)-qt0(i,j,k))/rdt
          e12p(i,j,k) = (1-fnudgeloc(i,j,k))*e12p(i,j,k) + fnudgeloc(i,j,k)*(e120save(i,j,k,1)-e120(i,j,k))/rdt
        end do
      end do
    end do

  end subroutine nudgeboundary

  subroutine exitnudgeboundary

    if (.not. (lnudgeboundary)) return
    deallocate(fnudgeglob,fnudgeloc)
    deallocate(umsave, vmsave, wmsave, thlmsave, qtmsave, e12msave)
    deallocate(u0save, v0save, w0save, thl0save, qt0save, e120save)
    deallocate(ql0save, ql0hsave, dthvdzsave,ekmsave,tmp0save,eslsave,qvslsave,qvsisave,presfsave,preshsave, thv0hsave)
    deallocate(u0avsave,v0avsave,thl0avsave,qt0avsave, thvhsave)
    deallocate(ksfcsave)
  end subroutine exitnudgeboundary

  function iglob(iloc,myidxloc)
    use modglobal, only : imax

    implicit none
    integer iloc,iglob,myidxloc

    iglob = iloc + imax*myidxloc - 1

  end function iglob

  function jglob(jloc,myidyloc)
    use modglobal, only : jmax

    implicit none
    integer jloc,jglob,myidyloc

    jglob = jloc + jmax*myidyloc - 1

  end function jglob

end module modnudgeboundary 

