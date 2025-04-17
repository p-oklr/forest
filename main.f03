program main
	use models									! Inflationary model
	use diagnostics
	use mpi
	use parameters								! Some parameters and math constants
	use random									! Custom random library
	implicit none
	
	!****************************************************!
	! ---------------- About PBHs -----------------------!
	!****************************************************!
	real :: pbh_threshold

	!****************************************************!
	! ---------------- About Stochasticity --------------!
	!****************************************************!

	! TODO Make use of the number of dimensions
	integer, parameter :: ndim = 1				! Number of dimensions
	real, parameter :: split_time = log(2.) / 3.0

	! Details of the sampling
	integer :: max_recursion					! Max level of recursion
	integer :: n_trees							! Number of independent trees
	integer :: max_tree							! Number of tree for current process

	! To reproduce individual trees
	logical :: replay
	integer, dimension(4) :: replay_seed=0

	!****************************************************!
	! --------------- Other variables -------------------!
	!****************************************************!
	class (abstract_model), allocatable :: model

	! MPI variables
	integer(4) :: rank, size_Of_Cluster, ierror, root

	!****************************************************!
	!-------------------- MAIN --------------------------!
	!****************************************************!

	! Initialize MPI
	call MPI_INIT(ierror)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, size_Of_Cluster, ierror)
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
	root = 0

	! Initialize the parametera of the simulation
	call init_program()

	! Initialize the diagnostics
	call init_diagnostics()

	! Manual load balancing
	max_tree = n_trees / size_Of_Cluster
	if (rank < MOD(n_trees, size_Of_Cluster)) then
		max_tree = max_tree + 1
	end if

	! Main loop
	call main_loop(max_tree)

	! Communication between the processes, reduction to the root node
	call reduce_diagnostics(rank, size_Of_Cluster)

	! Print the results stored in the root node
	if (rank == root) then 
		call print_results()
	end if

	! End MPI 
	call MPI_FINALIZE(ierror)

CONTAINS

!***************************************************************!
!---------- Initializes the program before execution -----------!
!***************************************************************!
subroutine init_program()

	integer(4) :: io								! Unit number, to read files
	character(len=64) :: model_type				! To receive the name of the inflation model

	! Recover parameter values at runtime using a namelist
	namelist /CONFIG/ max_recursion, n_trees, model_type, pbh_threshold, replay, replay_seed

	open (unit=io, file="config.nml")
	read (nml=CONFIG, unit=io)
	close (io)

	! If in replay mode, then only one tree
	if (replay) then 
		n_trees = 1
	end if

	! Initialize the parameters of the model
	if (model_type == "phi_square") then
		allocate (phi_square :: model)
	else if (model_type == "quantum_well") then
		allocate (quantum_well :: model)
	else if (model_type == "usr") then
		allocate (usr :: model)
	else
		print *, "Choose an existing model"
	end if

	call model%init_model()

end subroutine init_program

!***************************************************************!
!---------- Main do loop, parallelized using OpenMP ------------!
!***************************************************************!
subroutine main_loop(max_tree)

	integer, intent(in) :: max_tree
	integer :: i								! Dummy integer
	integer :: n								! Size of the random seed
	real :: volume								! Volume of the current tree
	real :: x_variable							! X-variable as in Eq. (A.1)

	integer, allocatable :: seed(:)				! Random seed of the PRNG
	real, allocatable :: pbhs(:)				! PBHs formed in the current tree

	! Allocate the random seed
	call random_seed(size=n)
	allocate(seed(n))

	do i=1, max_tree

		if (replay) then
			call random_seed(put=replay_seed)
		end if

		! Registers the seed of this given tree
		call random_seed(get=seed)
		call reset_stamp()

		! Explores the tree and recovers its physical properties
		call explore_tree(0, model%phi_0, volume, x_variable, pbhs, 1)

		! Registers its properties in a diagnostic
		call register(seed, volume, x_variable, pbhs)

	end do

end subroutine main_loop




!***************************************************************!
!-----  Implementation of the Euler-Maruyama's scheme ----------!
!***************************************************************!
real function langevin_euler_maruyama(field, delta_t, efold) result(phi_new)

	real, intent(in) :: field					! Value of the field
	real, intent(in) :: delta_t					! Time increment
	real, intent(in) :: efold					! Number of efolds

	! Updates the former field value with the drift and noise term
	! Over a time delta_t (in units of efolds)
	phi_new = field &
		- model%v_p_over_v(field, efold) * delta_t  &
		+ (sqrt(model%potential(field) * delta_t / 3.0) / (2.0 * pi)) * random_stdnormal()

	! Checks if it crosses the reflective boundary
	if (phi_new > model%phi_mirror) then
		phi_new = 2.0 * model%phi_mirror - phi_new
	endif

end function langevin_euler_maruyama


!***************************************************************!
!---------  Routine to explore the stochastic tree -------------!
!***************************************************************!
recursive subroutine explore_tree(nb_in, phi_in, volume, x_variable, pbhs, code)

	integer, intent(in) :: nb_in					! Level of recursion
	real, intent(in) :: phi_in						! Initial value of the field
	integer, intent(in) :: code

	real, intent(out) :: volume						! Volume of the current branch
	real, intent(out) :: x_variable					! X-variable of the current branch
	real, intent(out), allocatable :: pbhs(:)		! PBHs formed in this branch

	real :: phi_new									! Updated field value
	real :: efold									! Efold of the node (if it is a leaf)
	real :: gamma_1, gamma_2						! Amplitude between Laplacian and coarse-shelled
	real :: d2_zeta_1, d2_zeta_2					! Laplacian of the curvature perturbation

	real, dimension(:), allocatable :: pbhs_1, pbhs_2	! PBHs formed in lower branches
	real :: volume_1, volume_2							! Volume in lower branches 
	real :: x_variable_1, x_variable_2						! X-variable in lower branches

	real :: time									! Number of efolds on this node
	real :: delta_t									! Variable step size in efolds

	! At worse, return an allocated empty array
	pbhs = [real ::]

	!****************************************************!
	!-------- Too many levels of recursion --------------!
	!****************************************************!

	if ((nb_in >= max_recursion)) then
		stop "Too many levels of recursion"

	!****************************************************!
	!------------ Exploration of the tree ---------------!
	!****************************************************!

	else
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Propagate until the next splitting !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		time = 0.0
		phi_new = phi_in

		! Compute the current number of efolds
		efold = nb_in * split_time

		do while ((model%test_rho(phi_new) .eqv. .false.) .and. (time < split_time))

			! Get the best step size
			delta_t = model%get_delta_t(phi_new, min_efold, efold)

			! Do not go beyond the split_time
			delta_t = min(delta_t, split_time - time)

			! Update time
			time = time + delta_t
			efold = efold + delta_t

			! Update the value of the field
			phi_new = langevin_euler_maruyama(phi_new, delta_t, efold)			
		end do


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! If a crossing happened, register volume and efold number !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if (model%test_rho(phi_new) .eqv. .true.) then

			! Compute volume and X-variable
			volume = exp(3.0 * time)
			x_variable = volume * efold

			call register_leaf(volume, efold)

			if (replay) then
				print *, code, efold
			end if


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Otherwise, split the node into several branches !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		else


			! Explore the left branch
			call explore_tree(nb_in + 1, phi_new, volume_1, x_variable_1, pbhs_1, 2 * code)


			! Explore one of the right branch
			call explore_tree(nb_in + 1, phi_new, volume_2, x_variable_2, pbhs_2, 2 * code + 1)

			! Update local values
			volume = volume_1 + volume_2
			x_variable = x_variable_1 + x_variable_2

			!!!!!!!!!!!!!!!!!
			! Look for PBHs !
			!!!!!!!!!!!!!!!!!
			! Compute the coarse-shelled curvature perturbations
			gamma_1 = 2.0 / log(volume / volume_1)
			d2_zeta_1 = gamma_1 * (x_variable_1 / volume_1 - x_variable / volume)

			gamma_2 = 2.0 / log(volume / volume_2)
			d2_zeta_2 = gamma_2 * (x_variable_2 / volume_2 - x_variable / volume)

			! Check if it is above the threshold
			if (d2_zeta_1 > pbh_threshold) then 
				pbhs = [volume_1]

				if (replay) then
					print *, 2 * code, -1
				end if

			else
				pbhs = pbhs_1
			end if

			if (d2_zeta_2 > pbh_threshold) then 
				pbhs = [pbhs, [volume_2]]

				if (replay) then
					print *, 2 * code + 1, -1
				end if

			else
				pbhs = [pbhs, pbhs_2]
			end if
		endif
	endif

end subroutine explore_tree

end program main
