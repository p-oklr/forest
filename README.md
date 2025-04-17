# FOREST: FOrtran Recursive Exploration of Stochastic Trees

## How to cite

The code is made public on Github under the GPLv3 License. Its releases are referenced in Zenodo. 
If you use or modify the code, please cite the associated Zenodo DOI (either the version or concept DOI) as well as the accompanying publication:
```
@article{Animali:2025pyf,
    author = "Animali, Chiara and Auclair, Pierre and Blachier, Baptiste and Vennin, Vincent",
    title = "{Harvesting primordial black holes from stochastic trees with $\texttt{FOREST}$}",
    eprint = "2501.05371",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    month = "1",
    year = "2025"
}
```

## How to use

### Installation

The main code is written in standard **Fortran (>= 2003)**, parallelized using **MPI**, and the plots are generated using **Python3**.
The attached Makefile takes care of everything for you.

To compile and run the code, just type
```bash
make
```
Under the hood, the code will create repositories, compile all the `.f03` files, link them and produce the binary `main.bin`.

### Running the code

A configuration file sets the inflationary model and all the parameters of the simulation.
The configuration file is read at runtime, therefore you do not need to recompile the code every time.
In details, the configuration file `config.nml` contains NAMELISTs and looks like the following
```
&CONFIG
max_recursion=10000,
n_trees=1000,
model_type="quantum_well",
pbh_threshold=1.0
replay=.false.
replay_seed=0
/

&DIAGNOSTICS
n_bins=200,
max_volume=1e3,
min_efold=1e-3
max_efold=1e1,
max_pbh=5e2,
n_pbhs=10,
debug=.false.
/

&MODEL_QUANTUM_WELL
phi_0=1.0,
phi_mirror=1.0,
phi_star=0.0,
mu=0.8,
d=0.0
/

```
The results of the simulation will be stored has `.dat` files in the `data/` folder.


## Implementation details

### Project structure

Once you clone the project,you may find the following files:
- `diagnostics.f03`: takes care of storing/writing the relevant pieces of information from our trees
- `main.f03`: main part of the program
- `models.f03`: contains the different inflationary models
- `parameters.f03`: some mathematical constants
- `random.f03`: custom module to generate gaussian random numbers

The project is divided into a number of different directories:
- `data/`: contains the raw output materials after execution
- `figures/`: contains output figures
- `plot_scripts`: contains the Python3 plot scripts for tree reconstruction

After execution, you may obtain different types of files:
- `hist_*.dat`: datafiles containing histograms for interesting quantities
- `hist_*.pdf`: figures produced by the plot script
- `*.o *.mod *.bin`: byproducts of the compilation, you may erase them

To clean up after the execution, you may use
```bash
make clean
```

### Parallelization with MPI

The Fortran code is parallelizable by construction, we distribute our population of stochastic trees evenly to the different processes using MPI and use `MPI_REDUCE` to aggregate the results in the root process.
The different MPI processes are independent by default.

We tested the scaling of the code on `Betelgeuse` with $10^7$ samples.
The machine has $128$ physical cores and $256$ available threads.
The scaling is satisfactory.
| Number of cores | Execution time | Theoretical bound | Ratio |
|-----------------|----------------|-------------------|-------|
|   1 | 04:23 | 04:23 | 1.00 |
|   2 | 02:59 | 02:11 | 1.36 |
|   4 | 01:30 | 01:05 | 1.37 |
|   8 | 00:45 | 00:32 | 1.37 |
|  16 | 00:23 | 00:16 | 1.42 |
|  32 | 00:12 | 00:08 | 1.48 |
|  64 | 00:06 | 00:04 | 1.60 |
| 128 | 00:04 | 00:02 | 2.06 |
| 256 | 00:02 | 00:01 | 2.57 |

### DEBUG mode and REPLAY

After exploring a tree, all of its characteristics are sent for registration to the `diagnostics` module using the subroutine `register()`.
The purpose of this module is to store this information as efficiently as possible.
Nonetheless, we *cannot* store all this information and we design our diagnostics to keep only the most significant parts.

By activating the **DEBUG** mode in the configuration file, you force the code to write all the information to the standard output: the random seed leading to the tree, the volume, the $X$ variable and the list of all its PBHs.
```fortran
if (debug) then
	print *, "Seed:", seed
	print *, volume, x_variable, pbhs
end if
```
This gives you a direct access to the population of trees that you sampled
```bash
 Seed: -8108577098746041401  6732522176003023725  8639711273546409784  1706616359944309566
   62.993140323056430        210.52879071098701        61.822242175725272     
 Seed: -4234113764249552013  4567178560778395455  1548424405124481271  2547140088162218362
   1.1820786557689964        6.5910523497475310E-002
 Seed: -7740717987550496533  4306623956446807911 -5093175040809322841  -808786586031999030
   1.6700688807003021       0.28550655396906505     
 Seed: -3880012574634973134  7796300060026243019  2507977735375367313  5117952001553427977
   47.116431546506831        99.790938409300963        39.873850378393605        6.1235125902843528
```
For instance, here are four tree realizations, the first has one PBH while the fourth has two.

**NOTE**: The DEBUG mode slows the code down significantly and can rapidly lead to writing files of several GB. It should be used with care.


If there is a realization that caught your eye, you can use the **REPLAY** mode to reexplore it and store its inner structure.
The REPLAY mode does the following:
- Feeds the `replay_seed` in the configuration file to the pseudo-random number generator
- Overrides `n_trees`, only one tree will be explored
- Writes to the standard output the list of all the leaves with their number of efolds since the root node, as well as all the nodes leading to PBHs.

If you pipe the output of the REPLAY mode to a file, you can the `plot_scripts/reconstruct_*.py` Python files to visualize these trees.
```bash
./main.bin > queue.txt
python plot_scripts/reconstruct_tree.py
python plot_scripts/reconstruct_map.py
```

In case you do not have the required Python modules to produce the plots, you can install them easily using pip.
```bash
python3 -m venv env
source env/bin/activate
pip install -r requirements.txt
```


## Physical inputs

### Inflationary models

As of now, we have three different inflationary models
- $V(\phi) = M^4 \phi^2$, `phi_square`
- $V(\phi) = 24 \pi^2 v_0 (1 + d \phi)$, `quantum_well`
- Ultra Slow-Roll model, more on that later, `usr`

Each of these models is implemented as a class in the module of `models.f03`.
A model *extends* `abstract_model` (available in Fortran since 2003), an abstract class that specifies what a model should look like.
Using abstract classes has a number of advantages for our application:

- Some functions are identical for the different models. Now we only write them once, the new models will immediatly inherit these procedures.
- It allows us to change inflationary model at runtime in the config file, without having to recompile, since the structure of the module is strict.

As an example, model `phi_square` is set by declaring these **required** elements
```fortran
!*****************************************************
!----------------- MODEL PHI SQUARE -----------------!
!*****************************************************
type, extends(abstract_model) :: phi_square
contains
    procedure :: init_model => quantum_well_init_model
	procedure :: potential => phi_square_potential
	procedure :: v_p_over_v => phi_square_v_p_over_v
end type phi_square

(...)
contains
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
```

**NOTE**: We allow the drift term to depend on the absolute number of efold.
This may come as a surprise, since one expects the drift term $V'(\phi) / V(\phi)$ to depend only on the field value.
This *trick* allows us to use some models of Ultra Slow-Roll with a non-vanishing momentum $\bar{\pi}$, which decays exponentially with the number of efolds and acts as an additional drift term.
More details can be found in:

> *Ultra-slow-roll inflation with quantum diffusion*, Chris Pattison, Vincent Vennin, David Wands, Hooshyar Assadullahi, JCAP 04 (2021) 080

### PBH formation

#### Formation criterion 

The formation of primordial black holes usually takes place in regions of large curvature. 
However, using the value of the curvature perturbation as a criterion for PBH formation is not always well-justified, in particular when quantum diffusion plays an important role. 
The reason is that $\zeta$ is defined with respect to a given observer, and relative to a ``background'' that encloses their whole observable universe. 
However, gravitational collapse into a black hole is a local process, and should not depend on averages performed over regions much larger than the Hubble radius when it takes place. 
In other words, $\zeta_i$ is not a fully local quantity since it is affected by the volume produced in regions of the tree that are very distant from node $i$, via their contributions to $W_1$. Using $\zeta_i$ to decide whether or not the region emerging from node $i$ collapses into a black hole thus seems problematic. 

For this reason, other quantities are often considered for PBH criteria, such as the comoving density contrast, or its non-perturbative generalisation, the compaction function.
The reason why the comoving density contrast is better suited than the curvature perturbation is that, through Poisson equation, it is proportional to the Laplacian of the curvature perturbation, hence in Fourier space it features an additional $k^2$ factor. 
This suppresses large-scale contributions in the coarse-grained density contrast, which is thus a more ``local'' tracer than the curvature perturbation.

The same idea can be implemented in stochastic trees along the lines of the ``coarse-shelled'' curvature perturbation proposed in
> *Statistics of coarse-grained cosmological fields in stochastic inflation*, Yuichiro Tada, Vincent Vennin, JCAP 02 (2022) 02

Consider a given node $i$, giving birth to two daughter nodes $\ell$ and $m$. Our goal is to use the node $i$ as a local background for the node $\ell$, and to evaluate the curvature perturbation in node $\ell$ relative to its local environment $i$,
$$
\zeta_{\ell i} = \zeta_\ell-\zeta_i
$$
This can be related to the Laplacian of the curvature perturbation, or more generally to the compaction function.
Eventually, we find (see accompagnying paper)
$$
\frac{\zeta_{\ell i}}{\zeta_{\ell i, \mathrm{c}}} = \frac{2}{\ln(V_i/V_\ell)} \left(\mathcal{N}_{i\to \ell} + W_\ell - W_i \right)
$$
In its current state, the code is perfectly able to handle cloud-in-cloud.

#### Control parameters

The configuration files hosts three variables specific to the PBH mass function
- `pbh_max`: the maximum volume of the region that form PBHs that we want to track in our histograms
- `pbh_threshold`: multiplicative factor in front of the PBH formation threshold (should be unity)


