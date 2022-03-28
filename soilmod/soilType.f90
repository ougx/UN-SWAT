
#define soilname    this%name
#define itype       this%isoiltype
#define Ksat        this%parval(1)
#define thetaS      this%parval(2)
#define thetaR      this%parval(3)
#define alpha       this%parval(4)
#define soln        this%parval(5)
#define solm        this%parval(6)
#define hentry      this%parval(4)

#define nearSat(h)  abs((h) - merge(hentry, 0., itype>=soilRetentionFunction%BrooksCorey)) < veryTiny
#define nearOne(s)  abs((s) - 1.0) < veryTiny
#define sr_vg(h)    (1. + (-alpha * h) ** soln) ** (-solm)
#define kr_vg(sr)   sr ** 0.5 * (1. - (1. - sr ** (1./solm)) ** solm) ** 2
#define sr_bc(h)    (h / hentry) ** (-soln)
#define kr_bc(sr)   sr ** (2.5 + 2. / soln)
#define sr_gn(h)    exp(0.5 * alpha * h) * (1. - 0.5 * alpha * h) ** (2. / (soln + 2.))
#define kr_gn(h)    exp(alpha * h)


module mod_soiltype
  use constants, only         :  NA, veryTiny, soilRetentionFunction, dh

  type m_soiltype
    ! model - parameters
    ! vg: Ksat thetaS thetaR
    character(len=:),allocatable  :: name
    integer                       :: isoiltype
    real, allocatable             :: parval(:)
  contains
    procedure                     :: initialize               ! read parameter into parvals (give a string)
    procedure                     :: saturation2moisture      ! calculate saturation based on moisture
    procedure                     :: moisture2saturation      ! calculate moisture based on saturation
    procedure                     :: pressure2saturation      ! calculate saturation based on pressure
    procedure                     :: saturation2conductivity  ! calculate saturation based on pressure
    procedure                     :: pressure2conductivity    ! calculate hydraulic conductivity based on pressure
    procedure                     :: pressure2capacity        ! calculate hydraulic conductivity based on pressure
  end type


  type(m_soiltype),   allocatable :: soiltype(:)
contains

  subroutine initialize(this, line)
    class(m_soiltype)             :: this
    character(len=*)              :: line

    !local
    integer                       :: npar, iitype
    ! do nothing
    read(line, *) iitype
    select case (iitype)
    case (soilRetentionFunction%vanGenuchten)
      npar = 6
    case (soilRetentionFunction%BrooksCorey)
      npar = 6
    case (soilRetentionFunction%Gardner)
      npar = 6
    case default
      print*, 'Undefined Soil Model Type at Line:'//adjustl(trim(line))
      stop
    end select

    this%isoiltype = iitype
    allocate(this%parval(npar))
    read(line, *) iitype, this%parval
  end subroutine

  function saturation2moisture(this, saturation)
    class(m_soiltype)         :: this
    real                      :: saturation
    real                      :: saturation2moisture
    saturation2moisture = thetaR + (thetaS - thetaR) * saturation
  end function

  function moisture2saturation(this, moisture)
    class(m_soiltype)         :: this
    real                      :: moisture
    real                      :: moisture2saturation
    moisture2saturation = (moisture - thetaR) / (thetaS - thetaR)
  end function


  function pressure2saturation(this, pressure)
    class(m_soiltype)         :: this
    real                      :: pressure
    real                      :: pressure2saturation

    if (nearSat(pressure)) then
      pressure2saturation = 1.0
    else
      select case (itype)
      case (soilRetentionFunction%vanGenuchten)
        pressure2saturation = sr_vg(pressure)
      case (soilRetentionFunction%BrooksCorey)
        pressure2saturation = sr_bc(pressure)
      case (soilRetentionFunction%Gardner)
        pressure2saturation = sr_gn(pressure)
      end select
    end if

  end function

  function saturation2conductivity(this, saturation)
    class(m_soiltype)         :: this
    real                      :: saturation
    real                      :: saturation2conductivity

    ! local
    real                      :: kr

    if (nearOne(saturation)) then
      saturation2conductivity = Ksat
    else
      select case (itype)
      case (soilRetentionFunction%vanGenuchten)
        kr = kr_vg(saturation)
      case (soilRetentionFunction%BrooksCorey)
        kr = kr_bc(saturation)
      case (soilRetentionFunction%Gardner)
        print *, 'Kr is calculaed based on pressure in Gardner method '
        stop
      end select
      saturation2conductivity = Ksat * kr
    end if
  end function

  function pressure2capacity(this, pressure)
    class(m_soiltype)         :: this
    real                      :: pressure
    real                      :: pressure2capacity

    if (nearSat(pressure)) then
      pressure2capacity = 0.0
    else
      select case (itype)
      case (soilRetentionFunction%vanGenuchten)
        pressure2capacity = (thetaS - thetaR) * alpha * soln * solm * (1. + (- alpha * pressure) ** soln) ** (-solm-1.) * (- alpha * pressure) ** (soln - 1.)
      case (soilRetentionFunction%BrooksCorey)
        pressure2capacity = sr_bc(pressure)
      case (soilRetentionFunction%Gardner)
        pressure2capacity = sr_gn(pressure)
      end select
    end if

  end function


  function pressure2conductivity(this, pressure)
    class(m_soiltype)         :: this
    real                      :: pressure
    real                      :: pressure2conductivity

    ! local
    real                      :: sr, kr

    if (nearSat(pressure)) then
      pressure2conductivity = Ksat
    else
      select case (itype)
      case (soilRetentionFunction%vanGenuchten)
        sr = sr_vg(pressure)
        kr = kr_vg(sr)
      case (soilRetentionFunction%BrooksCorey)
        sr = sr_bc(pressure)
        kr = kr_bc(sr)
      case (soilRetentionFunction%Gardner)
        kr = kr_gn(pressure)
      end select
      pressure2conductivity = Ksat * kr
    end if
  end function

  function dKdh(this, pressure, k0)
    ! fintie difference of dK/dh
    class(m_soiltype)         :: this
    real                      :: pressure
    real                      :: k0
    real                      :: dKdh

    ! local
    real                      :: h1, k1

    if (nearSat(pressure)) then
      dKdh = 0.
    else
      h1 = pressure + dh
      k1 = this%pressure2conductivity(h1)
      dKdh = (k1 - k0) / dh
    end if
  end function


end module
