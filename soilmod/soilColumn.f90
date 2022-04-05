
#define nnz         this%nz
#define ddz         this%dz
#define h_new       this%pressure
#define ssink       this%sink
#define satr        this%saturation
#define conc        this%concentration
#define capa        this%capacity
#define cond        this%conductivity
#define stype(iz)   soiltype(this%isoiltype(iz))
#define scap(iz)    (stype(iz)%parval(2) - stype(iz)%parval(3))
#define disp(iz)    stype(iz)%parval(8)

! TODO todoKavg
#define todoKavg    soilInterblock_Kavg%SpaceWeighted
#define darcy(K, hup, hdn, dz) (K) - (K) * ((hdn) - (hup)) / (dz)
#define upstream(kup, kdn)  this%k_inter(iz, kup, kdn, soilInterblock_Kavg%Upstream)

#define satr0       this%satr_strt
#define conc0       this%conc_strt
#define wflow       this%w_flow
#define wext        this%w_ext
#define wsto        this%w_sto
#define sflow       this%s_flow
#define sext        this%s_ext
#define ssto        this%s_sto


module mod_soilcolumn
  use constants,    only : soilInterblock_Kavg, soilBC, dh, verySmall, Courant_max, minSatSink, minConcFrac, ZERO, ONE, TWO, HALF
  use mod_soiltype, only : soiltype
  use solver

  type m_soilcolumn
    integer                             :: nz                 ! number of computing nodes
    integer                             :: ncomp              ! number of chemical species
    real, allocatable                   :: dz(:)              ! grid size          [L]
    integer, allocatable                :: isoiltype(:)       ! soil type          [--]
    real,    allocatable                :: pressure(:)        ! pressure           [L]
    real,    allocatable                :: concentration(:,:) ! concentration      [M/L3]
    real,    allocatable                :: sink(:)            ! source and sink    [1/T]
    real,    allocatable                :: capacity(:)        ! moisture capacity  [1/L]
    real,    allocatable                :: saturation(:)      ! saturation         [--]
    real,    allocatable                :: satr_strt(:)       ! saturation at the start of simulation [L]
    real,    allocatable                :: conc_strt(:,:)     ! concentration at the start of simulation [M/L3]
    real,    allocatable                :: conductivity(:)    ! hydraulic cond.    [L/T]
    ! ------------ for mass balance check, stored as volumn or length
    real,    allocatable                :: w_flow(:)          ! water w_flow                      [L]
    !real,    allocatable                :: w_sto(:)           ! storage change                    [L]
    real,    allocatable                :: w_ext(:)           ! actual external source and sink,  [L]
    ! ------------ for mass balance check, stored as mass
    real,    allocatable                :: s_flow(:,:)        ! solute transport w_flow           [M]
    !real,    allocatable                :: s_sto(:,:)         ! storage change                    [M]
    real,    allocatable                :: s_ext(:,:)         ! actual external source and sink   [M]
  contains
    procedure                           :: initialize
    procedure                           :: update_state
    procedure                           :: distrub_cond
    procedure                           :: k_inter
    procedure                           :: solve => solve_newton
    procedure                           :: solve_picard
    procedure                           :: solve_solute
    procedure                           :: reset_budget
    procedure                           :: write_waterbudget
    procedure                           :: write_solutebudget
  end type

  type(m_soilcolumn), allocatable :: soilcolumn(:)
contains



  subroutine initialize(this, nz, dz, isoil, press_init, ncomp, conc_init)
    class(m_soilcolumn)                 :: this
    integer, intent(in)                 :: nz
    integer, intent(in)                 :: isoil(nz)
    real, intent(in)                    :: dz(nz)
    real, intent(in)                    :: press_init(nz)
    integer, intent(in), optional       :: ncomp
    real, intent(in), optional          :: conc_init(:,:)

    ! allocate and assign
    this%nz = nz
    allocate(ddz, source=dz)
    allocate(this%isoiltype, source=isoil)
    allocate(h_new, source=press_init)
    allocate(satr(nz))
    allocate(capa(nz))
    allocate(cond(nz))
    allocate(ssink(nz))

    !allocate(wsto(nz))
    allocate(wext(nz))
    allocate(wflow(0:nz))


    if (present(ncomp)) then
      this%ncomp = ncomp
      if (ncomp>0) then
        allocate(this%conc_strt,     source=conc_init)
        allocate(this%concentration, source=conc_init)
        !allocate(ssto(nz,ncomp))
        allocate(sext(nz,ncomp))
        allocate(sflow(0:nz,ncomp))
      end if
    else
      this%ncomp = 0
    end if
    ! update initial state and store the initial saturation
    call this%update_state()
    allocate(satr0, source=satr)

    call this%reset_budget()
  end subroutine

  subroutine reset_budget(this)
    class(m_soilcolumn)                 :: this

    wext  = ZERO
    wflow = ZERO
    !wsto  = ZERO
    satr0 = satr

    if (this%ncomp>0) then
      !ssto = ZERO
      sext = ZERO
      sflow = ZERO
      conc0 = conc
    end if

  end subroutine

  subroutine write_waterbudget(this, time, ifile)
    use iso_fortran_env, only : output_unit
    class(m_soilcolumn)                 :: this
    real   , intent(in)                 :: time
    integer, optional                   :: ifile

    ! local
    integer                             :: iifile, iz
    real                                :: residual, sink_tot, schg_tot, in_tot, out_tot, schg

    iifile = output_unit
    if (present(ifile)) iifile = ifile
    sink_tot = ZERO
    schg_tot = ZERO

    ! write header
    write(iifile, '(A)') ' Time     Layer       Pressure        FlowTop        FlowBot           Sink        Storage          Error'
    do iz =1, nnz
      sink_tot = sink_tot + wext(iz)
      schg = (satr(iz) - satr0(iz)) * scap(iz) * ddz(iz)
      schg_tot = schg_tot + schg
      residual = wflow(iz-1) - wflow(iz) + wext(iz) - schg
      write(iifile, '(F10.3,I5,100ES15.7)') time, iz, h_new(iz), wflow(iz-1), wflow(iz), wext(iz), -schg, residual
    end do

    in_tot = ZERO
    out_tot = ZERO
    if (wflow(0)>ZERO) then
      in_tot = in_tot + wflow(0)
    else
      out_tot = out_tot - wflow(0)
    end if
    if (wflow(nnz)<ZERO) then
      in_tot = in_tot - wflow(nnz)
    else
      out_tot = out_tot + wflow(nnz)
    end if
    in_tot   = in_tot  + sum(pack(wext, wext>ZERO))
    out_tot  = out_tot - sum(pack(wext, wext<ZERO))

    if (schg_tot>ZERO) then
      out_tot = out_tot + schg_tot
    else
      in_tot  = in_tot  - schg_tot
    end if

    residual = wflow(0) - wflow(nnz) + sink_tot - schg_tot
    write(iifile, '(A,100ES15.7)') '              Total:', wflow(0), wflow(nnz), sink_tot, -schg_tot, residual
    write(iifile, *)  'Total In                          :', in_tot
    write(iifile, *)  'Total Out                         :', out_tot
    write(iifile, *)  'Mass balance discrepancy (percent):', 100 * (in_tot - out_tot) / in_tot

  end subroutine

  subroutine write_solutebudget(this, time, icomp, ifile)
    use iso_fortran_env, only : output_unit
    class(m_soilcolumn)                 :: this
    real   , intent(in)                 :: time
    integer                             :: icomp
    integer, optional                   :: ifile

    ! local
    integer                             :: iifile, iz
    real                                :: residual, sink_tot, schg_tot, in_tot, out_tot, schg

    iifile = output_unit
    if (present(ifile)) iifile = ifile
    sink_tot = ZERO
    schg_tot = ZERO

    ! write header
    write(iifile, '(A)') ' Time     Layer  Concentration        FlowTop        FlowBot           Sink        Storage          Error'
    do iz =1, nnz
      sink_tot = sink_tot + sext(iz, icomp)
      schg = (stype(iz)%saturation2moisture(satr(iz)) * conc(iz, icomp) - stype(iz)%saturation2moisture(satr0(iz)) * conc0(iz, icomp)) * ddz(iz)
      schg_tot = schg_tot + schg
      residual = sflow(iz-1, icomp) - sflow(iz, icomp) + sext(iz, icomp) - schg
      write(iifile, '(F10.3,I5,100ES15.7)') time, iz, conc(iz, icomp), sflow(iz-1, icomp), sflow(iz, icomp), sext(iz, icomp), -schg, residual
    end do

    in_tot = ZERO
    out_tot = ZERO
    if (wflow(0)>ZERO) then
      in_tot = in_tot + sflow(0, icomp)
    else
      out_tot = out_tot - sflow(0, icomp)
    end if
    if (wflow(nnz)<ZERO) then
      in_tot = in_tot - sflow(nnz, icomp)
    else
      out_tot = out_tot + sflow(nnz, icomp)
    end if
    in_tot   = in_tot  + sum(pack(sext(:, icomp), sext(:, icomp)>ZERO))
    out_tot  = out_tot - sum(pack(sext(:, icomp), sext(:, icomp)<ZERO))

    if (schg_tot>ZERO) then
      out_tot = out_tot + schg_tot
    else
      in_tot  = in_tot  - schg_tot
    end if

    residual = sflow(0, icomp) - sflow(nnz, icomp) + sink_tot - schg_tot
    write(iifile, '(A,100ES15.7)') '              Total:', sflow(0, icomp), sflow(nnz, icomp), sink_tot, -schg_tot, residual
    write(iifile, *)  'Total In                          :', in_tot
    write(iifile, *)  'Total Out                         :', out_tot
    write(iifile, *)  'Mass balance discrepancy (percent):', 100 * (in_tot - out_tot) / in_tot

  end subroutine

  function k_inter(this, iz, kup, kdn, ikavg)   ! k mean between iz and iz + 1
    class(m_soilcolumn)                 :: this
    integer, intent(in)                 :: iz
    real, intent(in)                    :: kup, kdn
    integer, intent(in), optional       :: ikavg
    real                                :: k_inter

    ! local
    integer                             :: iikavg

    if (present(ikavg)) then
      iikavg = ikavg
    else
      iikavg = soilInterblock_Kavg%SpaceWeighted ! todo: now interblock K is fixed
    end if

    select case (iikavg)
    case (soilInterblock_Kavg%SpaceWeighted)
      k_inter = (ddz(iz) * kup + ddz(iz+1) * kdn) / (ddz(iz) + ddz(iz+1))
    case (soilInterblock_Kavg%Upstream)
      if (h_new(iz)+(ddz(iz) + ddz(iz+1))/2 > h_new(iz+1)) then
        k_inter = kup
      else
        k_inter = kdn
      end if
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

  subroutine solve_newton(this, dt_sum, itopBC, topBC_val, ibotBC, botBC_val, t0, topConc)
    ! F. Hassane Maina and P. Ackerer, 2017: Ross scheme, Newton¨CRaphson iterative methods
    class(m_soilcolumn)                 :: this

    integer, intent(in)                 :: itopBC, ibotBC
    real   , intent(in)                 :: dt_sum, topBC_val
    real   , intent(in), optional       :: botBC_val      ! not needed if it is a free drainage
    real   , intent(in), optional       :: t0             ! starting time stamp (default is 0)
    real   , intent(in), optional       :: topConc        ! influx concentration from top

    ! local
    logical                             :: converge
    integer                             :: iz, iter, iter_tot, iter_nl, icomp
    real                                :: s_old(nnz), h_old(nnz), dhdt(nnz), dhdt_old(nnz), k_dif(nnz) !, h_tmp(nnz)
    real                                :: dzflow(0:nnz)    ! flow length           [cm]
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
    real                                :: ttopConc
    !integer                             :: iCrit
    integer                             :: info
    ! results
    real                                :: aflow(0:nnz)     ! actual flow
    real                                :: asink(nnz)       ! actual sink
    real                                :: aschg(nnz)       ! atual storage change


    ddt = dt_init
    dt_factor = ONE
    dhdt = ZERO
    iter_tot = 0
    iter_nl = 0
    if (present(t0)) then
      ct = t0
      ddt_sum = t0 + dt_sum
    else
      ct = ZERO
      ddt_sum = dt_sum
    end if
    if (present(topConc)) then
      ttopConc = topConc
    else
      ttopConc = ZERO
    end if

    ! calculate the flow path
    dzflow(0) = ddz(1)
    dzflow(1:nnz-1) = (ddz(1:nnz-1) + ddz(2:nnz))
    dzflow(nnz) = ddz(nnz)
    dzflow = dzflow / TWO


    timeLoop: do while(ct < ddt_sum)
#ifdef debugging
      print*, 'soil Current Time', ct
#endif
      h_old = h_new
      s_old = satr
      dhdt_old = dhdt


      timeStepping: do
        ! check if it very close to the end
        ddt = min(ddt * dt_factor, dt_max)
        if (ddt_sum - (ct + ddt) < verySmall .or. ct + ddt > ddt_sum) ddt = ddt_sum - ct

        ! reduce the sink if cell is dry
        do iz = 1, nnz
          if (scap(iz) * (s_old(iz) - minSatSink) + ssink(iz) * ddt < 0.) then
            ! avoid sink when the cell is very dry
            asink(iz) = ssink(iz)
          else
            asink(iz) = ZERO
          end if
        end do

        nonLinear: do iter = 1, maxiter
          iter_nl = iter_nl + 1

          ! set boundary conditions
          call this%distrub_cond(k_dif)


          iz = 0
          select case (itopBC)
          case (soilBC%head)
            gradient = 1 + (topBC_val - h_new(1)) / dzflow(0)
            aflow(0) = cond(1) * gradient
            dkavgdhdn = (k_dif(1) - cond(1)) / dh
            dqdhup(0) = ZERO
            dqdhdn(0) = dkavgdhdn * gradient - cond(1) / dzflow(0)

          case (soilBC%flux)
            aflow(0) = topBC_val
            dqdhup(0) = ZERO
            dqdhdn(0) = ZERO
          end select

          do iz = 1, nnz - 1
            kavg = this%k_inter(iz, cond(iz), cond(iz+1))
            gradient = 1 + (h_new(iz) - h_new(iz+1)) / dzflow(iz)
            aflow(iz) = kavg * gradient

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
            aflow(iz) = kavg * gradient
            dkavgdhup = (k_dif(iz) - cond(iz)) / dh
            dqdhup(iz) = dkavgdhup * gradient + cond(iz) / dzflow(iz)
            dqdhdn(iz) = ZERO

          case (soilBC%flux)
            aflow(iz) = topBC_val
            dqdhup(iz) = ZERO
            dqdhdn(iz) = ZERO

          case (soilBC%free)
            aflow(iz) = cond(iz)
            dqdhup(iz) = (k_dif(iz) - cond(iz)) / dh
            dqdhdn(iz) = ZERO

          end select

          ! assemble matrix
          do iz=1, nnz
            ! right hand side
            aschg(iz) = ddz(iz) / ddt * scap(iz) * (satr(iz) - s_old(iz))
            B(iz) = (aflow(iz-1) - aflow(iz)) + ddz(iz) * asink(iz) - aschg(iz)
            ! diagonal
            DDH(iz) = dqdhup(iz) - dqdhdn(iz-1) + ddz(iz) / ddt * capa(iz)
            ! sub-diagonal for the upper node
            DL(iz) = -dqdhup(iz-1)
            ! super-diagonal for the lower node
            DU(iz) =  dqdhdn(iz)
          end do

          ! solve for dh
          converge = sum(B**2)<verySmall
          if (converge) then
            ! the current head is good enough
            dt_factor = rmax
            exit nonLinear ! jump to next step because system does not change
          else
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
              call this%write_waterbudget(ct)
              print *, 'Warning: DGTSV Failed with info =', info
#endif
              exit nonLinear ! ill-posed problem to solve.
            end if
          end if

        end do nonLinear

        ! check truncation error
        if (converge) then
          dhdt = (h_new - h_old) / ddt
          err_trunc_max = maxval(ddt / TWO * abs(dhdt - dhdt_old))
          tolerance = max(truncate_rel * maxval(abs(h_new)) + truncate_abs, EPS*100)
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
            dt_factor = HALF
          end if
          ! reset the pressure
          h_new = h_old
          call this%update_state()
        else
          print*, 'Failed to find non-linear solution at t = ', ct, ' dt', ddt, ' dtfactor', dt_factor
          stop
        end if
      end do timeStepping
      wflow = wflow + aflow * ddt
      wext  = wext + asink * ddt * ddz
      do icomp = 1, this%ncomp
        call this%solve_solute(icomp, ddt, ttopConc, s_old, aflow, asink, dzflow, ct)
      end do
      ct = ct + ddt
    end do timeLoop
  end subroutine

  subroutine solve_picard(this, dt_sum, itopBC, topBC_val, ibotBC, botBC_val, t0, topConc)
    ! Celia, M., Bouloutas, E., Zarba, R., 1990. WRR
    class(m_soilcolumn)                 :: this

    integer, intent(in)                 :: itopBC, ibotBC
    real   , intent(in)                 :: dt_sum, topBC_val
    real   , intent(in), optional       :: botBC_val      ! not needed if it is a free drainage
    real   , intent(in), optional       :: t0             ! starting time stamp (default is 0)
    real   , intent(in), optional       :: topConc        ! influx concentration from top

    ! local
    logical                             :: converge
    integer                             :: iz, iter, iter_tot, iter_nl, icomp
    real                                :: s_old(nnz), h_old(nnz), h_chg(nnz), dhdt(nnz), dhdt_old(nnz)
    real                                :: dzflow(0:nnz)    ! flow length           [cm]
    real                                :: ddt, ct, ddt_sum
    real                                :: kavg, kavg_dz, dzdt
    doubleprecision, dimension(nnz)     :: DDH, B, DL, DU
    doubleprecision                     :: dh_max
    !real                                :: err_trunc(nnz)
    real                                :: err_trunc_max
    real                                :: tolerance
    real                                :: dt_factor
    real                                :: ttopConc
    !integer                             :: iCrit
    integer                             :: info
    ! results
    real                                :: aflow(0:nnz)     ! actual flow
    real                                :: asink(nnz)       ! actual sink
    !real                                :: aschg(nnz)       ! atual storage change


    ddt = dt_init
    dt_factor = ONE
    dhdt = ZERO
    iter_tot = 0
    iter_nl = 0
    if (present(t0)) then
      ct = t0
      ddt_sum = t0 + dt_sum
    else
      ct = ZERO
      ddt_sum = dt_sum
    end if
    if (present(topConc)) then
      ttopConc = topConc
    else
      ttopConc = ZERO
    end if

    ! calculate the flow path
    dzflow(0) = ddz(1)
    dzflow(1:nnz-1) = (ddz(1:nnz-1) + ddz(2:nnz))
    dzflow(nnz) = ddz(nnz)
    dzflow = dzflow / TWO


    timeLoop: do while(ct < ddt_sum)
      h_old = h_new
      s_old = satr
      dhdt_old = dhdt
      converge = .false.


      timeStepping: do
        ! check if it very close to the end
        ddt = min(ddt * dt_factor, dt_max)
        if (ddt_sum - (ct + ddt) < verySmall .or. ct + ddt > ddt_sum) ddt = ddt_sum - ct
        ! reduce the sink if cell is dry
        do iz = 1, nnz
          if (scap(iz) * (s_old(iz) - minSatSink) + ssink(iz) * ddt < ZERO) then
            ! avoid sink when the cell is very dry
            asink(iz) = ssink(iz)
          else
            asink(iz) = ZERO
          end if
        end do

        nonLinear: do iter = 1, maxiter
          iter_nl = iter_nl + 1

          ! initialize coefficients
          DDH = ZERO; B = ZERO; DL = ZERO; DU = ZERO


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

            B(iz) = B(iz) + kavg + dzdt * (scap(iz) * (satr(iz) - s_old(iz)) - capa(iz) * h_new(iz)) - ddz(iz) * asink(iz)

            DDH(iz+1) = -kavg_dz
            DL(iz+1) = kavg_dz
            B(iz+1) = -kavg
          end do

          iz = nnz
          dzdt = ddz(iz) / ddt
          B(iz) = B(iz) + dzdt * (scap(iz) * (satr(iz) - s_old(iz)) - capa(iz) * h_new(iz)) - ddz(iz) * asink(iz)
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
            call this%write_waterbudget(ct)
            print *, 'Warning: DGTSV Failed with info =', info
#endif
            exit nonLinear ! ill-posed problem to solve.
          end if
        end do nonLinear

        ! check truncation error
        if (converge) then
          dhdt = (h_new - h_old) / ddt
          err_trunc_max = maxval(ddt / TWO * abs(dhdt - dhdt_old))
          tolerance = max(truncate_rel * maxval(abs(h_new)) + truncate_abs, EPS*100.)
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
            dt_factor = HALF
          end if
          ! reset the pressure
          h_new = h_old
          call this%update_state()
        else
          print*, 'Failed to find non-linear solution at t = ', ct
          stop
        end if
      end do timeStepping


      ! calculate flow at the end
      select case (itopBC)
      case (soilBC%head)
        aflow(0) = cond(1) * ((topBC_val - h_new(1))/dzflow(0) + ONE)
      case (soilBC%flux)
        aflow(0) = topBC_val
      end select
      do iz = 1, nnz - 1
        kavg = this%k_inter(iz, cond(iz), cond(iz+1))
        aflow(iz) = kavg * ((h_new(iz) - h_new(iz+1))/dzflow(iz) + ONE)
      end do
      iz = nnz
      select case (ibotBC)
      case (soilBC%head)
        aflow(iz) = cond(iz) * ((h_new(iz) - botBC_val)/dzflow(iz) + ONE)
      case (soilBC%flux)
        aflow(iz) = botBC_val
      case (soilBC%free)
        aflow(iz) = cond(iz)
      end select

      wflow = wflow + aflow * ddt
      wext  = wext + asink * ddt * ddz
      do icomp = 1, this%ncomp
        call this%solve_solute(icomp, ddt, ttopConc, s_old, aflow, asink, dzflow, ct)
      end do
      ct = ct + ddt
    end do timeLoop
  end subroutine

  subroutine solve_solute(this, icomp, dt_sum, topConc, s_old, aflow, asink, dzflow, t0)
    ! solve the advection dispersion equation after richards equation
    ! top boundary assume constant concentration during a infiltration event
    ! zero concentration for evaporation
    ! the method used is picard
    ! qup * cup - qdn cdn + scap(iz) * satr * dispersivity * mean(v1, v2)
    class(m_soilcolumn)                 :: this

    integer, intent(in)                 :: icomp          ! index of the species for the simulation (currently not supporting reaction)
    real   , intent(in)                 :: dt_sum         ! time step length (maybe subdivided depending on the corunt number)
    real   , intent(in)                 :: topConc        ! concentration of the applied water (irrigation or precipitation)
    real   , intent(in)                 :: dzflow(0:nnz)  ! flow length           [cm]
    real   , intent(in)                 :: aflow(0:nnz)   ! actual flow
    real   , intent(in)                 :: asink(nnz)     ! actual sink
    real   , intent(in)                 :: s_old(nnz)     ! old saturation
    real   , intent(in)                 :: t0             ! starting time stamp

    ! local
    !logical                             :: converge
    real                                :: Courant        ! Courant number = v * dt / ddz, this should smaller than 2
    real                                :: maxConc        ! Courant number = v * dt / ddz, this should smaller than 2
    real, parameter                     :: omega_disp = ONE ! dispersion time weight, c = omega_disp * c^{t} + (1- omega_disp) * c^{t+1}

    integer                             :: iz, iter, nstep
    real                                :: vel
    real                                :: topflux, topflow
    real                                :: botflow
    real                                :: qdisp(0:nnz)       ! dispersion flux
    real                                :: disp_dn(nnz)     ! disperssivity
    real                                :: ups(nnz)         ! upstream weight
    real                                :: csink(nnz)       ! contaminant sink (including lateral flow, plant root uptake (transpiration) etc)
    real                                :: ss_new(nnz)      ! saturation at the end of transport step
    real                                :: ss_old(nnz)      ! saturation at the end of transport step
    real                                :: c_old(nnz)       ! concentration at the beginning of transport step
    real                                :: c_chg(nnz)       ! concentration change between iterations
    real                                :: s_chg(nnz)       ! saturation change rate
    doubleprecision, dimension(nnz)     :: DDH, B, DL, DU, DDH0, B0, DL0, DU0, SDT
    real                                :: ct, ddt, cddt0, dc_max
    integer                             :: info

    ct = t0
    ddt = dt_sum

    ! calculate velocity and mositure change
    maxConc = max(maxval(conc(:,icomp)), topConc) * minConcFrac

    do iz = 1, nnz
      ss_old(iz) = stype(iz)%saturation2moisture(s_old(iz))
      ss_new(iz) = stype(iz)%saturation2moisture(satr(iz))
      ! check time step size
      if (conc(iz, icomp) > maxConc) then
        vel = (aflow(iz-1) + aflow(iz)) / (ss_old(iz) + ss_new(iz))
        Courant = abs(vel) * ddt / ddz(iz)
        if (Courant-Courant_max>verySmall) ddt = ddt * Courant_max / Courant
      end if
    end do

    nstep = ceiling(dt_sum / ddt)
    cddt0 = dt_sum / nstep
    ddt = cddt0

    s_chg = (ss_new - ss_old) / dt_sum

    ! calculate some variables that are unchanged over the transport step
    topflow = max(ZERO, aflow(0))
    botflow = max(ZERO, aflow(nnz))
    topflux = topConc * topflow
    csink = merge(asink, ZERO, asink<ZERO) * ddz ! only account for outflow (sinks)
    ups =  merge(ONE, ZERO, aflow(1:nnz)>ZERO)
    do iz = 1, nnz-1
      disp_dn(iz) = this%k_inter(iz, disp(iz), disp(iz+1))
    end do

    timeLoop: do while(ct < t0+dt_sum)
#ifdef debugging
      print*, ' transport Current Time', ct
      !print*, 'initial concentration', conc(1:nnz, icomp)
#endif
      c_old = conc(:, icomp)
      ! calculate the dispersion flux
      qdisp(0)       = ZERO !disp(1)* topflow * (topConc - c_old(1)) / dzflow(0)
      qdisp(1:nnz-1) = disp_dn(1:nnz-1) * aflow(1:nnz-1) * (c_old(1:nnz-1) - c_old(2:nnz)) / dzflow(1:nnz-1)

      DDH0 = ZERO; B0 = ZERO; DU0 = ZERO; DL0 = ZERO

      ! top boundary
      B0(1) = topflux + qdisp(0)
      do iz = 1, nnz-1
        B0(iz)   = B0(iz) - qdisp(iz) !  + ss_old(iz) * c_old(iz) * ddz(iz) / ddt
        DDH0(iz) = DDH0(iz) + s_chg(iz) * ddz(iz) + ups(iz) * aflow(iz) - csink(iz) ! saturation_old over dt
        DU0(iz) = (ONE - ups(iz)) * aflow(iz)

        B0(iz+1) = qdisp(iz)
        DL0(iz+1) = -ups(iz) * aflow(iz)
        DDH0(iz+1) = -(ONE - ups(iz)) * aflow(iz)
      end do

      ! bottom bounary
      DDH0(nnz) = DDH0(nnz) + botflow - csink(nnz)

      timeStepping: do
        nonLinear: do iter = 1, maxiter
          if (t0+dt_sum - (ct + ddt) < verySmall .or. ct-t0 + ddt > dt_sum) ddt = dt_sum+t0 - ct
          SDT = ss_old * ddz / ddt
          DL = DL0
          DDH = DDH0 + SDT
          DU = DU0
          B = B0 + SDT * c_old
          call dgtsv(nnz, 1, DL(2:nnz), DDH, DU, B, nnz, info)
          if (info == 0) then
            ! linear solution found
            c_chg = B - conc(:, icomp)
            conc(:, icomp) = B
            ! check max head change
            dc_max = maxval(abs(c_chg))
            if (dc_max<=cclose) then
              exit timeStepping ! if not converge loop next iteration
            else
#ifdef debugging
              print *, ' cclose failed transport at time ', ct, ' Iteration', iter, ' dc_max', dc_max, 'topflux', topflux
              !print *, 'DDH', DDH
              !print *, 'B', B
              !print *, 'DU', DU
              !print *, 'DL', DL
              !stop
#endif

            end if
          else
            exit nonLinear
          end if
        end do nonLinear

        ddt = HALF * ddt
        if (ddt<dt_limit) then
          print*, 'Minimum time step reach before transport is solved for species', icomp, ' at time ', ct, ' dt', ddt, ' dt0', cddt0
          print*, conc(:, icomp)
          stop
        end if

      end do timeStepping
      qdisp(0) = qdisp(0) + topflux
      qdisp(1:nnz-1) = qdisp(1:nnz-1) + aflow(1:nnz-1) * merge(conc(1:nnz-1, icomp), conc(2:nnz, icomp), aflow(1:nnz-1)>0)
      qdisp(nnz) = botflow * conc(nnz, icomp)  ! this is not dispersion flux, save for convinience when calculating budget

      sflow(:,icomp) = sflow(:,icomp) + qdisp * ddt
      sext (:,icomp) = sext (:,icomp) + csink * ddt

      ss_old = ss_old + s_chg * ddt
      ct = ct + ddt
      if (ddt<cddt0 .and. iter<10) then
        ddt = ddt * TWO
      end if
    end do timeLoop

  end subroutine
end module
