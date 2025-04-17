module random
	use parameters
	implicit none
	private

	public :: random_stdnormal
	public :: reset_stamp

	real, save :: stamp = 0.                       ! Used to optimize RNG

CONTAINS

	!*****************************************************!
	!-- Generates a uniformly distributed random number --!
	!*****************************************************!
	real function random_stduniform()
		real :: r

		call random_number(r)
		random_stduniform = 1 - r

	end function random_stduniform

	!*****************************************************!
	!----------- Resets the value of the stamp -----------!
	!*****************************************************!
	! This procedure is interesting when one wants to generate identical sequences of random numbers
	subroutine reset_stamp()
		stamp = 0.
	end subroutine reset_stamp

	!*****************************************************!
	!--- Generates a normally distributed random number --!
	!*****************************************************!
	! via Box-Muller transform
	real function random_stdnormal()
		real :: u1, u2, s

		if (stamp == 0.) then

			s = 2.
			do while (s >= 1.)
				u1 = 2 * random_stduniform() - 1
				u2 = 2 * random_stduniform() - 1
				s = u1 ** 2 + u2 ** 2
			end do

			random_stdnormal = u1 * sqrt(-2*log(s) / s)
			stamp = u2 * sqrt(-2*log(s) / s)

		else
			random_stdnormal = stamp
			stamp = 0.
		end if

	end function random_stdnormal
end module
