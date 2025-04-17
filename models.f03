module models
    use parameters
    implicit none
    private

    public :: abstract_model
    public :: phi_square
    public :: quantum_well
    public :: usr

    !******************************************************
    !------ Define an abstract class for this model ------!
    !******************************************************

    type, abstract :: abstract_model
        real :: phi_0
        real :: phi_mirror
        real :: phi_star
    contains
        procedure :: get_delta_t
        procedure :: test_rho
        procedure(init_model_shape), deferred :: init_model
        procedure(potential_shape), deferred :: potential
        procedure(v_p_over_v_shape), deferred :: v_p_over_v
    end type abstract_model

    ! We need to provide an interface, that is how the yet unspecified functions
    ! work and what they need
    interface
        ! Model initialization
        subroutine init_model_shape(this)
            import :: abstract_model
            class(abstract_model), intent(inout) :: this
        end subroutine init_model_shape

        ! Potential
        real function potential_shape(this, field)
            import :: abstract_model
            class(abstract_model), intent(in) :: this
            real, intent(in) ::  field
        end function potential_shape

        ! Computing (derivative_V / V)
        real function v_p_over_v_shape(this, field, efold)
            import :: abstract_model
            class(abstract_model), intent(in) :: this
            real, intent(in) :: field, efold
        end function v_p_over_v_shape

    end interface

    !*****************************************************
    !----------------- MODEL PHI SQUARE -----------------!
    !*****************************************************
    type, extends(abstract_model) :: phi_square
        real :: mass
    contains
        procedure :: init_model => phi_square_init_model
        procedure :: potential => phi_square_potential
        procedure :: v_p_over_v => phi_square_v_p_over_v
    end type phi_square

    !*****************************************************
    !---------------- MODEL QUANTUM WELL ----------------!
    !*****************************************************
    type, extends(abstract_model) :: quantum_well
        real :: mu
        real :: d
    contains
        procedure :: init_model => quantum_well_init_model
        procedure :: potential => quantum_well_potential
        procedure :: v_p_over_v => quantum_well_v_p_over_v
    end type quantum_well

    !*****************************************************
    !--------------------- MODEL USR --------------------!
    !*****************************************************
    type, extends(abstract_model) :: usr
        real :: mu
        real :: y_in
    contains
        procedure :: init_model => usr_init_model
        procedure :: potential => usr_potential
        procedure :: v_p_over_v => usr_v_p_over_v
    end type usr

CONTAINS

    !******************************************************
    !----------------- COMMON PROCEDURES -----------------!
    !******************************************************

    ! Find the best value for delta_t
    real function get_delta_t(this, field, min_efold, efold) result(delta_t)
        class(abstract_model), intent(inout) :: this
        real, intent(in) :: field
        real, intent(in) :: min_efold
        real, intent(in) :: efold

        ! 5 sigma probability to cross
        delta_t = min(min_efold, 3.0 * (2.0 * pi * (field - this%phi_star)) ** 2.0 / this%potential(field) / 5.0)

        ! Drift term bound
        if (this%v_p_over_v(field, efold) /= 0.0) then
            delta_t = min(delta_t, (field - this%phi_star) / this%v_p_over_v(field, efold) / 10.0)
        end if

    end function get_delta_t

    ! Test of the end of inflation
    logical function test_rho(this, field)
        class(abstract_model), intent(in) :: this
        real, intent(in) :: field

        if (field <= this%phi_star) then
            test_rho = .true.
        else
            test_rho = .false.
        end if

    end function test_rho

    !*****************************************************
    !----------------- MODEL PHI SQUARE -----------------!
    !*****************************************************

    ! Initialize using the config file
    subroutine phi_square_init_model(this)
        class(phi_square), intent(inout) :: this
        integer :: io
        real :: phi_0
        real :: phi_mirror
        real :: phi_star
        real :: mass

        ! Recover parameter values at runtime using a namelist
        namelist /MODEL_PHI_SQUARE/ phi_0, phi_mirror, phi_star, mass

        ! Open and read parameter values !
        open (unit=io, file="config.nml")
        read (nml=MODEL_PHI_SQUARE, unit=io)
        close (io)

        ! Feed the data structure
        this%phi_0 = phi_0
        this%phi_mirror = phi_mirror
        this%phi_star = phi_star
        this%mass = mass

    end subroutine phi_square_init_model

    ! Potential
    real function phi_square_potential(this, field) result(potential)
        class(phi_square), intent(in) :: this
        real, intent(in) ::  field
        potential = this%mass**4 * field**2

    end function phi_square_potential

    ! Computing (derivative_V / V)
    real function phi_square_v_p_over_v(this, field, efold) result(v_p_over_v)
        class(phi_square), intent(in) :: this
        real, intent(in) :: field, efold
        v_p_over_v = 2.0 / field

    end function phi_square_v_p_over_v

    !*****************************************************
    !---------------- MODEL QUANTUM WELL ----------------!
    !*****************************************************

    ! Initialize using the config file
    subroutine quantum_well_init_model(this)
        class(quantum_well), intent(inout) :: this
        integer :: io
        real :: phi_0
        real :: phi_mirror
        real :: phi_star
        real :: mu
        real :: d

        ! Recover parameter values at runtime using a namelist
        namelist /MODEL_QUANTUM_WELL/ phi_0, phi_mirror, phi_star, mu, d

        ! Open and read parameter values !
        open (unit=io, file="config.nml")
        read (nml=MODEL_QUANTUM_WELL, unit=io)
        close (io)

        ! Feed the data structure
        this%phi_0 = phi_0
        this%phi_mirror = phi_mirror
        this%phi_star = phi_star
        this%mu = mu
        this%d = d

    end subroutine quantum_well_init_model

    ! Definition of the potential
    real function quantum_well_potential(this, field) result(potential)
        class(quantum_well), intent(in) :: this
        real, intent(in) ::  field
        real :: alpha

        alpha = this%d * (this%phi_mirror - this%phi_star)
        potential = 24.0 * pi ** 2.0 * ((this%phi_mirror - this%phi_star) / this%mu) ** 2.0 * (1 + alpha * field)

    end function quantum_well_potential

    ! Computing (derivative_V / V)
    real function quantum_well_v_p_over_v(this, field, efold) result(v_p_over_v)
        class(quantum_well), intent(in) :: this
        real, intent(in) :: field, efold
        real :: alpha

        alpha = this%d * (this%phi_mirror - this%phi_star)
        v_p_over_v = alpha / (1 + alpha * field)

    end function quantum_well_v_p_over_v

    !*****************************************************
    !--------------------- MODEL USR --------------------!
    !*****************************************************

    ! Initialize using the config file
    subroutine usr_init_model(this)
        class(usr), intent(inout) :: this
        integer :: io
        real :: phi_0
        real :: phi_mirror
        real :: phi_star
        real :: mu
        real :: y_in

        ! Recover parameter values at runtime using a namelist
        namelist /MODEL_USR/ phi_0, phi_mirror, phi_star, mu, y_in

        ! Open and read parameter values !
        open (unit=io, file="config.nml")
        read (nml=MODEL_USR, unit=io)
        close (io)

        ! Feed the data structure
        this%phi_0 = phi_0
        this%phi_mirror = phi_mirror
        this%phi_star = phi_star

        this%y_in = y_in
        this%mu = mu

    end subroutine usr_init_model


    ! Definition of the potential
    real function usr_potential(this, field) result(potential)
        class(usr), intent(in) :: this
        real, intent(in) ::  field

        potential = 24.0 * pi ** 2.0 * ((this%phi_mirror - this%phi_star) / this%mu) ** 2.0

    end function usr_potential

    ! Computing (derivative_V / V)
    real function usr_v_p_over_v(this, field, efold) result(v_p_over_v)
        class(usr), intent(in) :: this
        real, intent(in) :: field, efold
        real :: pi_crit
        
        pi_crit = -3.0 * (this%phi_mirror - this%phi_star)

        v_p_over_v = -pi_crit * this%y_in * exp(-3.0 * efold)

    end function usr_v_p_over_v
end module
