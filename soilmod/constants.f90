module constants
  private
  real, public, parameter                 :: ZERO = 0.D0
  real, public, parameter                 :: HALF = 5.D-1
  real, public, parameter                 :: ONE  = 1.D0
  real, public, parameter                 :: TWO  = 2.D0
  real, public, parameter                 :: NA = -9e30
  real, public, parameter                 :: veryTiny  = 1e-20
  real, public, parameter                 :: verySmall = 1e-10
  real, public, parameter                 :: dh = 1e-6            ! increment for dK/dh calculation
  real, public, parameter                 :: Courant_max = 4.     ! maximum courant number
  real, public, parameter                 :: minConcFrac = 0.2    ! courant number only evaluated for cells with concentration larger than the minimaum concentration fraction times the maximum concentration
  real, public, parameter                 :: minSatSink  = 0.1    ! if the saturation of a layer is smaller than the sink then the sink is reduced to zero

  type m_soilBC
    integer                               :: head = 1             ! head boundary
    integer                               :: flux = 2             ! flux
    integer                               :: free = 3             ! free drainage (unit hydraulic gradient))
  end type

  type m_soilRetentionFunction
    integer                               :: vanGenuchten = 1     ! van Genuchten
    integer                               :: Gardner      = 3     ! Gardner
    integer                               :: BrooksCorey  = 7     ! BrooksCorey
  end type

  type m_soilInterblock_Kavg
    integer                               :: SpaceWeighted = 1    ! mean = (dx1 v1 + dx2 v2) / (dx1 + dx2)
    integer                               :: linearInterp  = 2    ! mean = (dx1 v2 + dx2 v1) / (dx1 + dx2)
    integer                               :: Arithmetic    = 3    ! mean = (v1 + v2) * 0.5
    integer                               :: Geometric     = 4    ! mean = sqrt(v1*v2)
    integer                               :: Harmonic      = 5    ! mean = (dx1 + dx2) / (dx1 / v1 + dx2 / v2)
    integer                               :: Upstream      = 6    ! mean = v1 if h1 > h2 else v2
    integer                               :: Baker         = 7    ! Baker, D.L., 2006. General validity of conductivity means in unsaturated flow models. J. Hydrol. Eng. 11 (6), 526¨C538.
  end type

  type(m_soilBC)               , public, parameter   :: soilBC = m_soilBC(1,2,3)
  type(m_soilRetentionFunction), public, parameter   :: soilRetentionFunction = m_soilRetentionFunction(1,3,7)
  type(m_soilInterblock_Kavg)  , public, parameter   :: soilInterblock_Kavg = m_soilInterblock_Kavg(1,2,3,4,5,6,7)
end module
