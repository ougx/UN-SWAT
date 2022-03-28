
#define nnz         this%nz
#define ddz         this%dz
#define h_new       this%pressure
#define ssink       this%sink
#define fflow       this%flow
#define satr        this%saturation
#define satr0       this%satr_strt
#define capa        this%capacity
#define cond        this%conductivity
#define eext        this%ext
#define stype(iz)   soiltype(this%isoiltype(iz))
#define scap(iz)    (stype(iz)%parval(2) - stype(iz)%parval(3))

! TODO todoKavg
#define todoKavg    soilInterblock_Kavg%SpaceWeighted
#define darcy(K, hup, hdn, dz) (K) - (K) * ((hdn) - (hup)) / (dz)


module mod_soilcolumn
  use constants,    only : soilInterblock_Kavg, soilBC, dh, verySmall
  use mod_soiltype, only : soiltype
  use solver

  type m_soilcolumn
    integer                             :: nz              ! number of computing nodes
    real, allocatable                   :: dz(:)           ! grid size          [L]
    integer, allocatable                :: isoiltype(:)    ! soil type          [--]
    real,    allocatable                :: pressure(:)     ! pressure           [L]
    real,    allocatable                :: sink(:)         ! source and sink    [1/T]
    real,    allocatable                :: capacity(:)     ! moisture capacity  [1/L]
    real,    allocatable                :: saturation(:)   ! saturation         [--]
    real,    allocatable                :: satr_strt(:)    ! saturation at the start of simulation [L]
    real,    allocatable                :: conductivity(:) ! hydraulic cond.    [L/T]
    ! ------------ for mass balance check, stored as volumn or length
    real,    allocatable                :: flow(:)         ! water flow             [L]
    real,    allocatable                :: sto(:)          ! storage change         [L]
    real,    allocatable                :: ext(:)          ! actual external source and sink,  [L]
  contains
    procedure                           :: initialize
    procedure                           :: update_state
    procedure                           :: distrub_cond
    procedure                           :: k_inter
    procedure                           :: solve
    procedure                           :: solve_picard
    procedure                           :: writebudget
  end type

  type(m_soilcolumn), allocatable :: soilcolumn(:)
contains

  subroutine initialize(this, nz, dz, isoil, press_init)
    class(m_soilcolumn)                 :: this
    integer, intent(in)                 :: nz
    integer, intent(in)                 :: isoil(nz)
    real, intent(in)                    :: dz(nz)
    real, intent(in)                    :: press_init(nz)

    ! allocate and assign
    this%nz = nz
    allocate(ddz, source=dz)
    allocate(this%isoiltype, source=isoil)
    allocate(h_new, source=press_init)
    allocate(satr(nz))
    allocate(capa(nz))
    allocate(cond(nz))
    allocate(ssink(nz))
    allocate(eext(nz))
    allocate(fflow(0:nz))

    eext = 0.
    fflow = 0.

    ! update initial state and store the initial saturation
    call this%update_state()
    allocate(satr0, source=satr)

#ifdef debugging
  print*, 'initialized column with', nz, 'layers'
#endif
  end subroutine

  subroutine writebudget(this, time, ifile)
    use iso_fortran_env, only : output_unit
    class(m_soilcolumn)                 :: this
    real   , intent(in)                 :: time
    integer, optional                   :: ifile

    ! local
    integer                             :: iifile
    real                                :: sink, residual, sink_tot, schg_tot, in_tot, out_tot, schg

    iifile = output_unit
    ddt = 1.
    if (present(ifile)) iifle = ifile
    sink_tot = 0.
    schg_tot = 0.

    ! write header
    write(iifile, '(A)') ' Time     Layer       Pressure        FlowTop        FlowBot           Sink        Storage          Error'
    do iz =1, nnz
      sink = ddz(iz) * eext(iz)
      sink_tot = sink_tot + sink
      schg = (satr(iz) - satr0(iz)) * scap(iz) * ddz(iz)
      schg_tot = schg_tot + schg
      residual = fflow(iz-1) - fflow(iz) + sink - schg
      write(iifile, '(F10.3,I5,100ES15.7)') time, iz, h_new(iz), fflow(iz-1), fflow(iz), sink, -schg, residual
    end do

    in_tot = 0.
    out_tot = 0.
    if (fflow(0)>0.) then
      in_tot = in_tot + fflow(0)
    else
      out_tot = out_tot - fflow(0)
    end if
    if (fflow(nnz)<0.) then
      in_tot = in_tot - fflow(nnz)
    else
      out_tot = out_tot + fflow(nnz)
    end if
    in_tot   = in_tot  + sum(pack(ssink * ddz, ssink>0.))
    out_tot  = out_tot - sum(pack(ssink * ddz, ssink<0.))

    if (schg_tot>0) then
      out_tot = out_tot + schg_tot
    else
      in_tot  = in_tot  - schg_tot
    end if

    residual = fflow(0) - fflow(nnz) + sink_tot - schg_tot
    write(iifile, '(A,100ES15.7)') '              Total:', fflow(0), fflow(nnz), sink_tot, -schg_tot, residual
    write(iifile, *)  'Total In                          :', in_tot
    write(iifile, *)  'Total Out                         :', out_tot
    write(iifile, *)  'Mass balance discrepancy (percent):', 100 * (in_tot - out_tot) / in_tot

  end subroutine

  function k_inter(this, iz, kup, kdn)   ! k mean between iz and iz + 1
    class(m_soilcolumn)                 :: this
    integer, intent(in)                 :: iz
    real, intent(in)                    :: kup, kdn
    real                                :: k_inter

    ! local
    integer                             :: ikavg = soilInterblock_Kavg%SpaceWeighted ! todo: now interblock K is fixed

    select case (ikavg)
    case (soilInterblock_Kavg%SpaceWeighted)
      k_inter = (ddz(iz) * kup + ddz(iz+1) * kdn) / (ddz(iz) + ddz(iz+1))
    end select
  end function

  subroutine update_state(this)
    class(m_soilcolumn)                 :: this

    ! local
    integer                             :: iz

    do iz = 1, nnz
      satr(iz) = stype(iz)%pressure2saturation(h_new(iz))
      cond(iz) = stype(iz)%saturation2conductivity(satr(iz))
      capa(iz) = stype(iz)%pressure2capacity(h_new(iz))
    end do
  end subroutine

  subroutine distrub_cond(this, Ktmp)
    class(m_soilcolumn)                 :: this
    real, dimension(:), intent(out)     :: Ktmp
    ! local
    integer                             :: iz

    do iz = 1, nnz
      Ktmp(iz) = stype(iz)%pressure2conductivity(h_new(iz)+dh)
    end do
  end subroutine

  subroutine solve(this, dt_sum, itopBC, topBC_val, ibotBC, botBC_val, t0)
    ! F. Hassane Maina and P. Ackerer, 2017: Ross scheme, Newton¨CRaphson iterative methods
    class(m_soilcolumn)                 :: this

    integer, intent(in)                 :: itopBC, ibotBC
    real   , intent(in)                 :: dt_sum, topBC_val
    real   , intent(in), optional       :: botBC_val      ! not needed if it is a free drainage
    real   , intent(in), optional       :: t0             ! starting time stamp (default is 0)

    ! local
    logical                             :: converge
    integer                             :: iz, iter, iter_tot, iter_nl
    real                                :: s_old(nnz), h_old(nnz), dhdt(nnz), dhdt_old(nnz), k_dif(nnz) !, h_tmp(nnz)
    real                                :: dzflow(0:nnz)    ! flow length           [cm]
    real                                :: flow(0:nnz)      ! flow
    real                                :: schg(nnz)        ! storage change
    real                                :: dqdhup(0:nnz)    ! differential flow w.r.t. head upper
    real                                :: dqdhdn(0:nnz)    ! differential flow w.r.t. head lower
    real                                :: ddt, ct, ddt_sum
    real                                :: kavg, dkavgdhup, dkavgdhdn, gradient
    doubleprecision, dimension(nnz)     :: DDH, B, DL, DU
    doubleprecision                     :: dh_max
    !real                                :: err_trunc(nnz)
    real                                :: err_trunc_max
    real                                :: tolerance
    real                                :: dt_factor
    !integer                             :: iCrit
    integer                             :: info


    ddt = dt_init
    dt_factor = 1.
    dhdt = 0.
    iter_tot = 0
    iter_nl = 0
    if (present(t0)) then
      ct = t0
      ddt_sum = t0 + dt_sum
    else
      ct = 0.
      ddt_sum = dt_sum
    end if

    ! calculate the flow path
    dzflow(0) = ddz(1)
    dzflow(1:nnz-1) = (ddz(1:nnz-1) + ddz(2:nnz))
    dzflow(nnz) = ddz(nnz)
    dzflow = dzflow / 2.


    timeLoop: do while(ct < ddt_sum)
#ifdef debugging
      print*, 'Current Time', ct
#endif
      h_old = h_new
      s_old = satr
      dhdt_old = dhdt
      converge = .false.


      timeStepping: do
        ! check if it very close to the end
        ddt = min(ddt * dt_factor, dt_max)
        if (ddt_sum - (ct + ddt) < verySmall .or. ct + ddt > ddt_sum) ddt = ddt_sum - ct

        nonLinear: do iter = 1, maxiter
          iter_nl = iter_nl + 1
          print*, 'Current Time', ct, ' Iter', iter

          ! set boundary conditions
          call this%distrub_cond(k_dif)

          iz = 0
          select case (itopBC)
          case (soilBC%head)
            gradient = 1 + (topBC_val - h_new(1)) / dzflow(0)
            flow(0) = cond(1) * gradient
            dkavgdhdn = (k_dif(1) - cond(1)) / dh
            dqdhup(0) = 0.
            dqdhdn(0) = dkavgdhdn * gradient - cond(1) / dzflow(0)

          case (soilBC%flux)
            flow(0) = topBC_val
            dqdhup(0) = 0.
            dqdhdn(0) = 0.
          end select

          do iz = 1, nnz - 1
            kavg = this%k_inter(iz, cond(iz), cond(iz+1))
            gradient = 1 + (h_new(iz) - h_new(iz+1)) / dzflow(iz)
            flow(iz) = kavg * gradient

            dkavgdhup = (this%k_inter(iz, k_dif(iz), cond(iz+1)) - kavg) / dh
            dkavgdhdn = (this%k_inter(iz, cond(iz), k_dif(iz+1)) - kavg) / dh

            dqdhup(iz) = dkavgdhup * gradient + kavg / dzflow(iz)
            dqdhdn(iz) = dkavgdhdn * gradient - kavg / dzflow(iz)
          end do

          iz = nnz
          select case (ibotBC)
          case (soilBC%head)
            kavg = cond(iz)
            gradient = 1 + (h_new(iz) - botBC_val) / dzflow(iz)
            flow(iz) = kavg * gradient
            dkavgdhup = (k_dif(iz) - cond(iz)) / dh
            dqdhup(iz) = dkavgdhup * gradient + cond(iz) / dzflow(iz)
            dqdhdn(iz) = 0.

          case (soilBC%flux)
            flow(iz) = topBC_val
            dqdhup(iz) = 0.
            dqdhdn(iz) = 0.

          case (soilBC%free)
            flow(iz) = cond(iz)
            dqdhup(iz) = (k_dif(iz) - cond(iz)) / dh
            dqdhdn(iz) = 0.

          end select

          ! assemble matrix
          do iz=1, nnz
            ! right hand side
            schg(iz) = ddz(iz) / ddt * scap(iz) * (satr(iz) - s_old(iz))
            B(iz) = (flow(iz-1) - flow(iz)) + ddz(iz) * ssink(iz) - schg(iz)
            ! diagonal
            DDH(iz) = dqdhup(iz) - dqdhdn(iz-1) + ddz(iz) / ddt * capa(iz)
            ! sub-diagonal for the upper node
            DL(iz) = -dqdhup(iz-1)
            ! super-diagonal for the lower node
            DU(iz) =  dqdhdn(iz)
          end do

          ! solve for dh
          call dgtsv(nnz, 1, DL(2:nnz), DDH, DU, B, nnz, info)

          if (info == 0) then ! linear solution found
            h_new = h_new + B
            call this%update_state()
            ! check max head change
            dh_max = maxval(abs(B))
            converge = (dh_max<=hclose)
            if (converge) then
              exit nonLinear ! if not converge loop next iteration
            else
#ifdef debugging
              print '(A,F10.1,A,I5,A,F10.1,A,F10.3,A,ES10.3)', ' Time ', ct, ' Iteration', iter, ' dt', ddt, ' dtfactor', dt_factor, ' dh_max', dh_max
#endif
            end if
          else
#ifdef debugging
            print '(2A,F10.1,A,I5,A,F10.1,A,F10.3)', new_line('a'), ' Time ', ct, ' Iteration', iter, ' dt', ddt, ' dtfactor', dt_factor
            call this%writebudget(ct)
            print *, 'Warning: DGTSV Failed with info =', info
#endif
            exit nonLinear ! ill-posed problem to solve.
          end if
        end do nonLinear

        ! check truncation error
        if (converge) then
          dhdt = B / ddt
          err_trunc_max = maxval(ddt / 2. * abs(dhdt - dhdt_old))
          tolerance = truncate_rel * maxval(abs(h_new)) + truncate_abs
          if (err_trunc_max <= tolerance) then
            ! converge and increase time step
            if (time_stepping) dt_factor = min(safety_factor * sqrt(tolerance/max(err_trunc_max, EPS)), rmax)
              !dt_factor = safety_factor * sqrt(tolerance/max(err_trunc_max, EPS))
            exit timeStepping
          end if
        end if
        if (time_stepping .and. ddt > dt_limit) then
          ! not converge, try to reduce time step if time stepping is enabled
          if (converge) then
            ! satisfy dh_max but not err_trunc_max
            dt_factor = max(safety_factor * sqrt(tolerance/max(err_trunc_max, EPS)), rmin)
          else
            ! not satisfy dh_max nor err_trunc_max
            dt_factor = 0.5
          end if
          ! reset the pressure
          h_new = h_old
          call this%update_state()
        else
          print*, 'Failed to find non-linear solution at t = ', ct
          stop
        end if
      end do timeStepping
      ct = ct + ddt
      fflow = fflow + flow * ddt
      eext  = eext + ssink * ddt
    end do timeLoop
#ifdef debugging
      print '(2A,F10.1,A,I5,A,F10.1,A,F10.3)', new_line('a'), ' Time ', ct, ' Iter_tot', iter_nl, ' dt', ddt, ' dtfactor', dt_factor
    call this%writebudget(ct)
#endif
  end subroutine

  subroutine solve_picard(this, dt_sum, itopBC, topBC_val, ibotBC, botBC_val, t0)
    ! Celia, M., Bouloutas, E., Zarba, R., 1990. WRR
    class(m_soilcolumn)                 :: this

    integer, intent(in)                 :: itopBC, ibotBC
    real   , intent(in)                 :: dt_sum, topBC_val
    real   , intent(in), optional       :: botBC_val      ! not needed if it is a free drainage
    real   , intent(in), optional       :: t0             ! starting time stamp (default is 0)

    ! local
    logical                             :: converge
    integer                             :: iz, iter, iter_tot, iter_nl
    real                                :: s_old(nnz), h_old(nnz), h_chg(nnz), dhdt(nnz), dhdt_old(nnz)
    real                                :: flow(0:nnz)      ! flow
    real                                :: schg(nnz)        ! storage change
    real                                :: dzflow(0:nnz)    ! flow length           [cm]
    real                                :: ddt, ct, ddt_sum
    real                                :: kavg, kavg_dz, dzdt
    doubleprecision, dimension(nnz)     :: DDH, B, DL, DU
    doubleprecision                     :: dh_max
    !real                                :: err_trunc(nnz)
    real                                :: err_trunc_max
    real                                :: tolerance
    real                                :: dt_factor
    !integer                             :: iCrit
    integer                             :: info


    ddt = dt_init
    dt_factor = 1.
    dhdt = 0.
    iter_tot = 0
    iter_nl = 0
    if (present(t0)) then
      ct = t0
      ddt_sum = t0 + dt_sum
    else
      ct = 0.
      ddt_sum = dt_sum
    end if

    ! calculate the flow path
    dzflow(0) = ddz(1)
    dzflow(1:nnz-1) = (ddz(1:nnz-1) + ddz(2:nnz))
    dzflow(nnz) = ddz(nnz)
    dzflow = dzflow / 2.


    timeLoop: do while(ct < ddt_sum)
#ifdef debugging
      print*, 'Current Time', ct
#endif
      h_old = h_new
      s_old = satr
      dhdt_old = dhdt
      converge = .false.


      timeStepping: do
        ! check if it very close to the end
        ddt = min(ddt * dt_factor, dt_max)
        if (ddt_sum - (ct + ddt) < verySmall .or. ct + ddt > ddt_sum) ddt = ddt_sum - ct

        nonLinear: do iter = 1, maxiter
          iter_nl = iter_nl + 1
          print*, 'Current Time', ct, ' Iter', iter

          ! initialize coefficients
          DDH = 0.; B = 0.; DL = 0.; DU = 0.

          ! set boundary conditions
          select case (itopBC)
          case (soilBC%head)
            kavg_dz = cond(1)/dzflow(0)
            B(1) = -cond(1)-kavg_dz*topBC_val
            DDH(1) =  -kavg_dz
          case (soilBC%flux)
            B(1) = -topBC_val
          end select

          do iz = 1, nnz - 1
            kavg = this%k_inter(iz, cond(iz), cond(iz+1))
            kavg_dz = kavg / dzflow(iz)
            dzdt = ddz(iz) / ddt

            DDH(iz) = DDH(iz) - kavg_dz - dzdt * capa(iz)
            DU(iz) = kavg_dz

            B(iz) = B(iz) + kavg + dzdt * (scap(iz) * (satr(iz) - s_old(iz)) - capa(iz) * h_new(iz)) - ddz(iz) * ssink(iz)

            DDH(iz+1) = -kavg_dz
            DL(iz+1) = kavg_dz
            B(iz+1) = -kavg
          end do

          iz = nnz
          dzdt = ddz(iz) / ddt
          B(iz) = B(iz) + dzdt * (scap(iz) * (satr(iz) - s_old(iz)) - capa(iz) * h_new(iz)) - ddz(iz) * ssink(iz)
          DDH(iz) = DDH(iz) - dzdt * capa(iz)
          select case (ibotBC)
          case (soilBC%head)
            kavg_dz = cond(iz)/dzflow(iz)
            B(iz) = B(iz) - kavg_dz * botBC_val + cond(iz)
            DDH(iz) = DDH(iz) - kavg_dz

          case (soilBC%flux)
            B(iz) = B(iz) + botBC_val

          case (soilBC%free)
            B(iz) = B(iz) + cond(iz)

          end select

          ! solve for dh
          call dgtsv(nnz, 1, DL(2:nnz), DDH, DU, B, nnz, info)

          if (info == 0) then ! linear solution found
            h_chg = B - h_new
            h_new = B
            call this%update_state()
            ! check max head change
            dh_max = maxval(abs(h_chg))
            converge = (dh_max<=hclose)
            if (converge) then
              exit nonLinear ! if not converge loop next iteration
            else
#ifdef debugging
              print '(A,F10.1,A,I5,A,F10.1,A,F10.3,A,ES10.3)', ' Time ', ct, ' Iteration', iter, ' dt', ddt, ' dtfactor', dt_factor, ' dh_max', dh_max
#endif
            end if
          else
#ifdef debugging
            print '(2A,F10.1,A,I5,A,F10.1,A,F10.3)', new_line('a'), ' Time ', ct, ' Iteration', iter, ' dt', ddt, ' dtfactor', dt_factor
            call this%writebudget(ct)
            print *, 'Warning: DGTSV Failed with info =', info
#endif
            exit nonLinear ! ill-posed problem to solve.
          end if
        end do nonLinear

        ! check truncation error
        if (converge) then
          dhdt = h_chg / ddt
          err_trunc_max = maxval(ddt / 2. * abs(dhdt - dhdt_old))
          tolerance = truncate_rel * maxval(abs(h_new)) + truncate_abs
          if (err_trunc_max <= tolerance) then
            ! converge and increase time step
            if (time_stepping) dt_factor = min(safety_factor * sqrt(tolerance/max(err_trunc_max, EPS)), rmax)
              !dt_factor = safety_factor * sqrt(tolerance/max(err_trunc_max, EPS))
            exit timeStepping
          end if
        end if
        if (time_stepping .and. ddt > dt_limit) then
          ! not converge, try to reduce time step if time stepping is enabled
          if (converge) then
            ! satisfy dh_max but not err_trunc_max
            dt_factor = max(safety_factor * sqrt(tolerance/max(err_trunc_max, EPS)), rmin)
          else
            ! not satisfy dh_max nor err_trunc_max
            dt_factor = 0.5
          end if
          ! reset the pressure
          h_new = h_old
          call this%update_state()
        else
          print*, 'Failed to find non-linear solution at t = ', ct
          stop
        end if
      end do timeStepping

      ct = ct + ddt

      ! calculate flow at the end
      select case (itopBC)
      case (soilBC%head)
        flow(0) = cond(1) * ((topBC_val - h_new(1))/dzflow(0) + 1)
      case (soilBC%flux)
        flow(0) = topBC_val
      end select
      do iz = 1, nnz - 1
        kavg = this%k_inter(iz, cond(iz), cond(iz+1))
        flow(iz) = kavg * ((h_new(iz) - h_new(iz+1))/dzflow(iz) + 1)
      end do
      iz = nnz
      select case (ibotBC)
      case (soilBC%head)
        flow(iz) = cond(iz) * ((h_new(iz) - botBC_val)/dzflow(iz) + 1)
      case (soilBC%flux)
        flow(iz) = botBC_val
      case (soilBC%free)
        flow(iz) = cond(iz)
      end select

      fflow = fflow + flow * ddt
      eext  = eext + ssink * ddt
    end do timeLoop
#ifdef debugging
      print '(2A,F10.1,A,I5,A,F10.1,A,F10.3)', new_line('a'), ' Time ', ct, ' Iter_tot', iter_nl, ' dt', ddt, ' dtfactor', dt_factor
      call this%writebudget(ct)
#endif
  end subroutine

end module
