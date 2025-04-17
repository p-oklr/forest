module diagnostics
    use mpi
    implicit none
    private

    public :: init_diagnostics
    public :: register
    public :: register_leaf
    public :: reduce_diagnostics
    public :: print_results

    public :: min_efold

    ! Parameters of the histogram
    integer(4) :: n_bins
    integer(4) :: n_pbhs
    integer :: n_trees = 0

    logical :: debug

    ! Tables for histograms, and their number of bins
    integer, dimension(:), allocatable :: vol_hist
	integer, dimension(:), allocatable :: efold_hist
	integer, dimension(:), allocatable :: pbh_hist
	integer, dimension(:, :), allocatable :: cluster_hist
	real, dimension(:), allocatable :: zeta_hist

    real :: total_volume = 0.0

    ! Maximum range for the histograms
	real :: max_volume							! Largest volume in histograms
	real :: min_efold							! Smallest efold in histograms
	real :: max_efold							! Largest efold in histograms
	real :: max_pbh								! Largest PBH mass

	! Bin widths for the histograms
	real :: log_bin_vol
	real :: log_bin_efold
	real :: log_bin_pbh

CONTAINS

!***************************************************************!
!----------------- Initializes the histograms ------------------!
!***************************************************************!
subroutine init_diagnostics()
    integer(4) :: io

    ! Recover the parameters needed to establish our diagnostics
    namelist /DIAGNOSTICS/ n_bins, n_pbhs, min_efold, max_efold, max_volume, max_pbh, debug
    open (unit=io, file="config.nml")
	read (nml=DIAGNOSTICS, unit=io)
	close (io)
	
    ! Allocate the tables for histograms
	allocate(vol_hist(n_bins))
	allocate(efold_hist(n_bins))
	allocate(pbh_hist(n_bins))
	allocate(cluster_hist(n_pbhs, n_bins))
	allocate(zeta_hist(n_bins))

    ! Zero the newly allocated tables
	vol_hist = 0.0
	efold_hist = 0.0
	pbh_hist = 0.0
	cluster_hist = 0.0
	zeta_hist = 0.0

    ! Initialize all variables based on the new configuration
	log_bin_vol = log(max_volume) / n_bins
	log_bin_efold = (log(max_efold) - log(min_efold)) / n_bins
	log_bin_pbh = log(max_pbh) / n_bins


end subroutine init_diagnostics

!***************************************************************!
!--------- Registers physical values into histograms -----------!
!***************************************************************!
subroutine register(seed, volume, x_variable, pbhs)

    integer, intent(in) :: seed(:)				! Random seed of the PRNG
	real, intent(in) :: volume					! Volume of the current tree
	real, intent(in) :: x_variable				! X-variable as in Eq. (A.1)
	real, intent(in) :: pbhs(:)					! PBHs formed in the current tree

	integer :: location							! Position in the histogram
	integer :: i								! Dummy integer

    ! Update global quantities
    total_volume = total_volume + volume
    n_trees = n_trees + 1

    if (debug) then
        print *, "Seed:", seed
        print *, volume, x_variable, pbhs
    end if

	!****************************************************!
	!------------------- Volume -------------------------!
	!****************************************************!

	! Bin number for the volume
	location = 1 + int(log(volume) / log_bin_vol)

	if ((location > 0) .and. (location < n_bins)) then
		vol_hist(location) = vol_hist(location) + 1
	end if

	!****************************************************!
	!------ Volume-averaged number of efolds ------------!
	!****************************************************!

	! Bin number for the average number of efolds
	location = 1 + int((log(x_variable / volume) - log(min_efold)) / log_bin_efold)

	if ((location > 0) .and. (location < n_bins)) then
		efold_hist(location) = efold_hist(location) + 1
	end if

	!****************************************************!
	!-------------------- PBHs --------------------------!
	!****************************************************!

	if (size(pbhs) > 0) then

		! Iterate over all the PBHs of the tree
		do i=1, size(pbhs)

			! Bin number for the PBH mass
			location = 1 + int(log(pbhs(i)) / log_bin_pbh)

			if ((location > 0) .and. (location < n_bins)) then
				pbh_hist(location) = pbh_hist(location) + 1
			end if
		enddo

	endif

	!****************************************************!
	!------ Number of PBHs per initial patch ------------!
	!****************************************************!

	if (size(pbhs) < n_pbhs) then 		
		! Bin number for the volume
		location = 1 + int(log(volume) / log_bin_vol)

		if ((location > 0) .and. (location < n_bins)) then
			cluster_hist(1 + size(pbhs), location) = cluster_hist(1 + size(pbhs), location) + 1
		end if
	end if

end subroutine

!***************************************************************!
!---- Registers physical values of leaves into histograms ------!
!***************************************************************!
subroutine register_leaf(volume, efold)

	real, intent(in) :: volume					! Volume of the current tree
	real, intent(in) :: efold				! X-variable as in Eq. (A.1)

	integer :: location							! Position in the histogram

	!****************************************************!
	!-------------- Number of efolds --------------------!
	!****************************************************!

	! Bin number for the average number of efolds
	location = 1 + int((log(efold) - log(min_efold)) / log_bin_efold)

	if ((location > 0) .and. (location < n_bins)) then
		zeta_hist(location) = zeta_hist(location) + volume
	end if


end subroutine


!***************************************************************!
!------------------ MPI Reduction to root node -----------------!
!***************************************************************!
subroutine reduce_diagnostics(rank, size_Of_Cluster)

    integer(4), intent(in) :: rank
    integer(4), intent(in) :: size_Of_Cluster

    ! MPI variables
    integer(4) :: count, ierror, root

    root = 0
    count = 1

    ! MPI Reduction to root node and print
	if (rank == root) then

		! Integer histograms
		call MPI_REDUCE(MPI_IN_PLACE, vol_hist, n_bins, MPI_INTEGER8, MPI_SUM, root, MPI_COMM_WORLD, ierror)
		call MPI_REDUCE(MPI_IN_PLACE, efold_hist, n_bins, MPI_INTEGER8, MPI_SUM, root, MPI_COMM_WORLD, ierror)
		call MPI_REDUCE(MPI_IN_PLACE, pbh_hist, n_bins, MPI_INTEGER8, MPI_SUM, root, MPI_COMM_WORLD, ierror)

		! 2D histogram
		call MPI_REDUCE(MPI_IN_PLACE, cluster_hist, n_bins * n_pbhs, MPI_INTEGER8, MPI_SUM, root, MPI_COMM_WORLD, ierror)

		! Real histogram
		call MPI_REDUCE(MPI_IN_PLACE, zeta_hist, n_bins, MPI_REAL8, MPI_SUM, root, MPI_COMM_WORLD, ierror)

        ! Integer number of trees in histogram
		call MPI_REDUCE(MPI_IN_PLACE, n_trees, count, MPI_INTEGER8, MPI_SUM, root, MPI_COMM_WORLD, ierror)

		! Real Volume
		call MPI_REDUCE(MPI_IN_PLACE, total_volume, count, MPI_REAL8, MPI_SUM, root, MPI_COMM_WORLD, ierror)
	else

		! Integer histogram
		call MPI_REDUCE(vol_hist, vol_hist, n_bins, MPI_INTEGER8, MPI_SUM, root, MPI_COMM_WORLD, ierror)
		call MPI_REDUCE(efold_hist, efold_hist, n_bins, MPI_INTEGER8, MPI_SUM, root, MPI_COMM_WORLD, ierror)
		call MPI_REDUCE(pbh_hist, pbh_hist, n_bins, MPI_INTEGER8, MPI_SUM, root, MPI_COMM_WORLD, ierror)

		! 2D histogram
		call MPI_REDUCE(cluster_hist, cluster_hist, n_bins * n_pbhs, MPI_INTEGER8, MPI_SUM, root, MPI_COMM_WORLD, ierror)

		! Real histogram
		call MPI_REDUCE(zeta_hist, zeta_hist, n_bins, MPI_REAL8, MPI_SUM, root, MPI_COMM_WORLD, ierror)

        ! Integer number of trees in histogram
		call MPI_REDUCE(n_trees, n_trees, count, MPI_INTEGER8, MPI_SUM, root, MPI_COMM_WORLD, ierror)

		! Real volume
		call MPI_REDUCE(total_volume, total_volume, count, MPI_REAL8, MPI_SUM, root, MPI_COMM_WORLD, ierror)
	endif

end subroutine reduce_diagnostics

!***************************************************************!
!------ Prints the results of the simulation to datafiles ------!
!***************************************************************!
subroutine print_results()

	integer :: i								! Dummy integer
	integer(4) :: io								! Unit number, to write files

	real :: bin_edge

	! Histogram for the volume
	open (unit=io, file="data/hist_vol.dat")
		write(io,*) 1, 0
		do i=1, n_bins
			bin_edge = exp(i * log_bin_vol)
			write(io, *) bin_edge, real(vol_hist(i)) / log_bin_vol / n_trees / bin_edge
		end do
	close (unit=io)

	! Histogram for the number of efolds
	open (unit=io, file="data/hist_zeta.dat")
		write(io,*) 0, 0
		do i=1, n_bins
			bin_edge = min_efold * exp(i * log_bin_efold)
			write(io, *) bin_edge, zeta_hist(i) / log_bin_efold / total_volume / bin_edge
		end do
	close (unit=io)


	! Histogram for the volume-averaged number of efolds
	open (unit=io, file="data/hist_efold.dat")
		write(io,*) 0, 0
		do i=1, n_bins
			bin_edge = min_efold * exp(i * log_bin_efold)
			write(io, *) bin_edge, real(efold_hist(i)) / log_bin_efold / n_trees / bin_edge
		end do
	close (unit=io)

	! Histogram for PBH formation
    open (unit=io, file="data/hist_pbhs.dat")
        write(io,*) 0, 0
        do i=1, n_bins
            bin_edge = exp(i * log_bin_pbh)
            write(io, *) bin_edge, real(pbh_hist(i)) / log_bin_pbh / total_volume / bin_edge
        end do
    close (unit=io)

    ! Histogram for the PBH clustering
    open (unit=io, file="data/hist_cluster.dat")
        do i=1, n_bins
            bin_edge = exp(i * log_bin_vol)
            write(io, *) bin_edge, real(cluster_hist(:, i)) / log_bin_vol / n_trees / bin_edge
        end do
    close (unit=io)

end subroutine print_results

end module diagnostics