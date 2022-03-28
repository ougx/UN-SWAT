module constants
  private
  real, public, parameter                 :: NA = -9e30
  real, public, parameter                 :: veryTiny  = 1e-20
  real, public, parameter                 :: verySmall = 1e-10
  real, public, parameter                 :: dh = 1e-6            ! increment for dK/dh calculation

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
    integer                               :: SpaceWeighted = 1    !
    integer                               :: Arithmetic    = 2    !
    integer                               :: Geometric     = 3    !
    integer                               :: Harmonic      = 4    !
    integer                               :: Upstream      = 5    !
    integer                               :: Baker         = 6    ! Baker, D.L., 2006. General validity of conductivity means in unsaturated flow models. J. Hydrol. Eng. 11 (6), 526¨C538.
  end type

  type(m_soilBC)               , public, parameter   :: soilBC = m_soilBC(1,2,3)
  type(m_soilRetentionFunction), public, parameter   :: soilRetentionFunction = m_soilRetentionFunction(1,3,7)
  type(m_soilInterblock_Kavg)  , public, parameter   :: soilInterblock_Kavg = m_soilInterblock_Kavg(1,2,3,4,5,6)
end module
