! Clonal simulation of the murine spermatogenic stem cell pool

! To compile the simulation, see compile.sh for examples how to invoke
! gfortran (required).

program germline
  implicit none

  ! constants and hard-wired parameters
  real(kind=8), parameter :: pi2 = 6.283185307179586
  real(kind=8), parameter :: h = 0.

  real(kind=8), parameter :: p_z2_labelling = 0.45

  character(255), parameter :: sim_directory = './'

  integer, parameter :: max_cs = 32
  integer, parameter :: max_n_div = 250

  integer, parameter :: num_layers = 20, num_syncytia = 8

  integer, parameter :: ALL = 0
  integer, parameter :: UNITS = 1, CELLS = 2

  ! functions
  integer :: rmod, unitstep, free_layer, occ_layer, seed

  ! parameters
  real(kind=8) :: nu, prob_diff_1, frag_1, div_2, frag_2, prob_dediff, conv_2, div_3, frag_3, diff_3, dt, sim_time, eq_time
  integer :: L_mig, L1, L2, nr, m

  ! derived parameters
  integer :: num_rps, eq_rp, num_sites, num_syncytia_2, num_ids
  real(kind=8) :: phi_g, diff_2, dediff_2, dediff_3

  ! fields
  integer, dimension(:, :), allocatable :: lt1, id1
  integer, dimension(:, :, :), allocatable :: lt2, lt3, id2, id3
  integer, dimension(:, :), allocatable :: pm1
  integer, dimension(:, :, :), allocatable :: pm2, pm3
  integer :: osc

  ! static arrays
  real(kind=8), dimension(:), allocatable :: equal_one, equal_two, larger_zero, larger_one, up_to_half

  ! observables
  real(kind=8) :: last_avg_lost_ids
  integer, dimension(:, :, :, :), allocatable :: csd, last_csd
  real(kind=8), dimension(:, :, :, :, :), allocatable   :: avg_csd
  real(kind=8), dimension(:), allocatable :: marg_avg_csd, avg_total_n_div_exit, avg_total_div_events
  real(kind=8), dimension(:, :, :), allocatable :: curdist
  real(kind=8), dimension(1:3) :: avg_content
  real(kind=8), dimension(0:3) :: origin, last_origin
  real(kind=8), dimension(1:3, 1:4) :: cell_flux, last_cell_flux
  integer, dimension(:), allocatable :: total_dupl_fates

  real(kind=8), dimension(:), allocatable :: lu, last_lu, lc, last_lc_frac, avg_lost_ids, avg_osc
  real(kind=8), dimension(:, :), allocatable :: avg_lu, avg_lc_frac, avg_n_div, avg_origin, avg_div_total
  integer, dimension(:, :), allocatable :: cmp, last_cmp
  real(kind=8), dimension(:), allocatable :: avg_lost_cells, last_div_total
  real(kind=8), dimension(:, :, :), allocatable :: avg_cmp, avg_cell_flux, avg_labelled_cmp
  real(kind=8), dimension(1:3) :: last_n_div
  real(kind=8), dimension(:, :), allocatable :: last_n_div_dist
  real(kind=8), dimension(:, :, :), allocatable :: avg_n_div_dist

  real(kind=8) :: lf2, last_lf2, lf3, last_lf3
  real(kind=8), dimension(:), allocatable :: times, avg_lf2, avg_lf3

  integer, dimension(:, :), allocatable :: num_units, last_num_units, num_marked_units, last_num_marked_units, num_labelled_units, last_num_labelled_units
  real(kind=8), dimension(:, :, :), allocatable :: avg_marked_ratio
  integer, dimension(:, :, :), allocatable :: nrm_marked_ratio

  ! computation variables
  real(kind=8), dimension(:), allocatable :: osc_signal
  integer :: unit_type

  real(kind=8), dimension(:, :, :), allocatable :: prop1
  real(kind=8), dimension(:, :, :, :), allocatable :: prop2
  real(kind=8), dimension(:, :, :, :), allocatable :: prop3
  real(kind=8) :: prop_osc, total_prop, t, time_step, reaction, p_sum, die, last_osc
  integer :: last_bp, orig_length, temp_lt, temp_id, inv
  logical :: searching, equilibrated
  integer :: temp_pm, lost_cells, last_lost_cells, num_cells, total_div_events, last_total_div_events, total_n_div_exit, last_total_n_div_exit

  ! counters
  integer :: i, a, layer, s, k, r, n, x1, x2, y1, y2, c, uc, rp, lrp, rpi, n1, n2, nt

  ! file names
  character(255) :: directory, directory_r, filename, filename_frame, title

  ! command line arguments
  integer, parameter :: req_num_args = 21
  integer :: num_args
  character(255), dimension(1:req_num_args) :: cl_args

  real(kind=8) :: start_time

  ! read command line arguments
  num_args = iargc()
  if (num_args.ne.req_num_args) then
    write(*,'(a)') 'Wrong number of arguments. Terminated.'
    call exit()
  end if
  do i = 1, num_args
    call getarg(i, cl_args(i))
  end do

  read(cl_args(1),*) nu  ! global oscillator transition rate
  read(cl_args(2),*) m  ! global oscillator number of states

  read(cl_args(3),*) div_2  ! incomplete division rate of the Y compartment
  read(cl_args(4),*) div_3  ! incomplete division rate of the Z compartment
  read(cl_args(5),*) dediff_3  ! Z -> Y transition rate

  read(cl_args(6),*) prob_diff_1  ! probability of X -> Y differentiation upon X division

  read(cl_args(7),*) frag_1  ! fragmentation rate of the X compartment
  read(cl_args(8),*) frag_2  ! fragmentation rate of the Y compartment

  read(cl_args(9),*) prob_dediff  ! probability of a Y -> X transition when Y unit converts
  read(cl_args(10),*) conv_2  ! conversion rate of the Y compartment

  read(cl_args(11),*) frag_3  ! fragmentation rate of the Z compartment
  read(cl_args(12),*) diff_3  ! differentiation rate of the Z compartment

  read(cl_args(13),*) L_mig  ! migration radius

  read(cl_args(14),*) L1  ! spatial extension in dimension 1
  read(cl_args(15),*) L2  ! spatial extension in dimension 2

  read(cl_args(16),*) dt  ! output time step

  read(cl_args(17),*) sim_time  ! total simulation time
  read(cl_args(18),*) eq_time  ! equilibration time

  read(cl_args(19),*) nr  ! number of model realisations

  read(cl_args(20),*) seed  ! RNG seed

  title = trim(cl_args(21))  ! name of the output directory

  ! create directory
  directory = trim(sim_directory) // trim(title) // '/'
  call system('mkdir -p ' // trim(directory))

  ! derived parameters
  num_rps = nint(sim_time / dt)
  eq_rp = nint(eq_time / dt)
  num_sites = L1 * L2
  num_ids = num_sites
  num_syncytia_2 = floor(real(num_syncytia, 8) / 2.)

  dediff_2 = prob_dediff * conv_2
  diff_2 = (1. - prob_dediff) * conv_2

  ! allocate dynamic arrays
  allocate(lt1(1:L1, 1:L2), lt2(1:num_layers, 1:L1, 1:L2), lt3(1:num_layers, 1:L1, 1:L2))
  allocate(id1(1:L1, 1:L2), id2(1:num_layers, 1:L1, 1:L2), id3(1:num_layers, 1:L1, 1:L2))
  allocate(pm1(1:L1, 1:L2), pm2(1:num_layers, 1:L1, 1:L2), pm3(1:num_layers, 1:L1, 1:L2))

  allocate(prop1(1:2, 1:L1, 1:L2), prop2(1:4, 1:num_layers, 1:L1, 1:L2), prop3(1:4, 1:num_layers, 1:L1, 1:L2))

  ! allocate observables
  allocate(csd(0:max_cs, 0:max_cs, 0:max_cs, 1:2), last_csd(0:max_cs, 0:max_cs, 0:max_cs, 1:2))
  allocate(avg_csd(0:num_rps, 0:max_cs, 0:max_cs, 0:max_cs, 1:2), marg_avg_csd(0:max_cs))
  allocate(curdist(0:max_cs, 0:max_cs, 0:max_cs))
  allocate(cmp(1:3, 0:num_syncytia), last_cmp(1:3, 0:num_syncytia), avg_cmp(0:num_rps, 1:3, 0:num_syncytia), avg_lost_ids(0:num_rps))
  allocate(avg_lf2(0:num_rps), avg_lf3(0:num_rps), lu(1:3), last_lu(1:3), avg_lu(0:num_rps, 1:3), lc(1:3), last_lc_frac(1:3), avg_lc_frac(0:num_rps, 1:3))
  allocate(num_units(1:3, 1:num_syncytia), num_marked_units(1:3, 1:num_syncytia), avg_marked_ratio(0:num_rps, 1:3, 1:num_syncytia), nrm_marked_ratio(0:num_rps, 1:3, 1:num_syncytia), num_labelled_units(1:3, 1:num_syncytia), last_num_labelled_units(1:3, 1:num_syncytia))
  allocate(last_num_units(1:3, 1:num_syncytia), last_num_marked_units(1:3, 1:num_syncytia), last_div_total(1:3))
  allocate(avg_lost_cells(0:num_rps), avg_n_div(0:num_rps, 1:3), avg_total_n_div_exit(0:num_rps), avg_labelled_cmp(0:num_rps, 1:3, 1:num_syncytia), avg_div_total(0:num_rps, 1:3), avg_total_div_events(0:num_rps))
  allocate(avg_n_div_dist(0:num_rps, 1:3, 0:max_n_div), last_n_div_dist(1:3, 0:max_n_div))
  allocate(avg_origin(0:num_rps, 0:3), avg_cell_flux(0:num_rps, 1:3, 1:4))
  allocate(osc_signal(0:m), avg_osc(0:num_rps), total_dupl_fates(1:2))

  allocate(equal_one(0:num_syncytia), equal_two(0:num_syncytia), larger_zero(0:num_syncytia), larger_one(0:num_syncytia), up_to_half(0:num_syncytia))

  ! initalize static arrays
  equal_one = 0.
  equal_one(1) = 1.
  equal_two = 0.
  equal_two(2) = 1.
  larger_zero = 1.
  larger_zero(0) = 0.
  larger_one = 1.
  larger_one(0:1) = 0.
  up_to_half = 0.
  up_to_half(1:num_syncytia_2) = 1.

  ! initialise random number generator
  call init_random_seed(seed)

  ! initialise global observables
  avg_osc = 0
  avg_csd = 0.
  avg_cmp = 0.
  avg_labelled_cmp = 0.
  avg_lf2 = 0.
  avg_lf3 = 0.
  avg_marked_ratio = 0.
  nrm_marked_ratio = 0
  avg_lu = 0.
  avg_lost_cells = 0
  avg_n_div = 0.
  avg_lost_ids = 0.
  avg_total_div_events = 0.
  avg_total_n_div_exit = 0.
  avg_n_div_dist = 0.
  avg_div_total = 0.
  avg_origin = 0
  avg_cell_flux = 0
  total_dupl_fates = 0

  ! initialise oscillator
  ! -- constant propensity
  prop_osc = real(m,8) * nu

  ! -- oscillatory signal
  do i = 0, m
    osc_signal(i) = (sqrt(0.5 * pi2) * gamma(1. + 0.5 * h) / gamma(0.5 + 0.5 * h)) * sin(0.5 * pi2 * real(i,8) / real(m,8)) ** h
  end do

  ! iterate over realizations
  do r = 1, nr
    write(*,'(a,f6.2,a)',advance='no') '\b\b\b\b\b\b\b', 100. * real(r - 1, 8) / real(nr, 8), '%'

    ! initialise timers and flags
    equilibrated = .false.
    t = -eq_time
    rp = -1
    lrp = -1

    ! physical initial condition
    lt1 = 1; lt2 = 0; lt3 = 0
    id1 = 0; id2 = 0; id3 = 0
    pm1 = 0; pm2 = 0; pm3 = 0
    osc = 0.

    ! time evolution
    do while (rp < num_rps)
      ! initialise clone IDs when the equilibration time is reached
      if ((t >= 0.).and.(equilibrated.eqv..false.)) then
        ! reset timer to zero
        rp = 0

        ! initialise clone IDs - label X layer only
        a = 0
        id2 = 0
        id3 = 0
        do x1 = 1, L1
          do x2 = 1, L2
#ifdef labelx
              ! label all x cells
            a = a + 1
            id1(x1, x2) = a
#endif

#ifdef labely
            ! label all y cells
            do layer = 1, num_layers
              if (lt2(layer, x1, x2) > 0) then
                a = a + 1
                id2(layer, x1, x2) = a
              end if
            end do
#endif

#ifdef labelz
            ! label all z cells
            do layer = 1, num_layers
              if (lt3(layer, x1, x2) > 0) then
                a = a + 1
                id3(layer, x1, x2) = a
              end if
            end do
#endif

#ifdef labelz2
            ! label fraction of z units
            do layer = 1, num_layers
              if (lt3(layer, x1, x2) > 0) then
                call random_number(die)
                if (die < p_z2_labelling) then
                  a = a + 1
                  id3(layer, x1, x2) = a
                end if
              end if
            end do
#endif

#ifdef labelcomp
            ! label all compartments with a compartment-specific label
            id1(x1, x2) = 1
            do layer = 1, num_layers
              if (lt2(layer, x1, x2) > 0) then
                id2(layer, x1, x2) = 2
              end if
              if (lt3(layer, x1, x2) > 0) then
                id3(layer, x1, x2) = 3
              end if
            end do
            a = 3
#endif
          end do
        end do
        num_ids = a

        ! reset timer
        t = 0.
        equilibrated = .true.

        ! - initialise proliferation markers
        pm1 = 0; pm2 = 0; pm3 = 0

        ! - counters
        lost_cells = 0
        total_div_events = 0
        total_n_div_exit = 0
        cell_flux = 0

        ! - determine clone sizes and cmps for each barcode
        call analyze_clones(max_cs, num_ids, num_layers, num_syncytia, L1, L2, lt1, lt2, lt3, id1, id2, id3, csd, cmp)

          ! - determine origins of cells
#ifdef labelcomp
        origin = 0
        do x1 = 1, L1
          do x2 = 1, L2
            origin(id1(x1,x2)) = origin(id1(x1,x2)) + lt1(x1,x2)
            do layer=1,num_layers
              origin(id2(layer, x1, x2)) = origin(id2(layer, x1, x2)) + lt2(layer, x1, x2)
              origin(id3(layer, x1, x2)) = origin(id3(layer, x1, x2)) + lt3(layer, x1, x2)
            end do
          end do
        end do
#endif

        ! - set first observations
        last_csd = csd
        last_cmp = cmp

        last_lf2 = 0
        last_lf3 = 0
        last_num_units = -1
        last_num_marked_units = 0
        last_num_labelled_units = 0
        last_lu = 0
        last_lc_frac = 0
        last_lost_cells = 0
        last_total_div_events = 0
        last_total_n_div_exit = 0
        last_div_total = 0
        last_origin = origin
      end if

      ! --- dynamics ---

      ! -- propensities

      ! compute partial propensities
      prop1 = 0.
      prop2 = 0.
      prop3 = 0.
      do x1 = 1, L1
        do x2 = 1, L2
          ! compartment 1 - process 1: incomplete division (triggered by oscillator)
          prop1(1, x1, x2) = 0.

          ! compartment 1 - process 2: fragmentation
          prop1(2, x1, x2) = larger_one(lt1(x1, x2)) * (lt1(x1, x2) - 1) * frag_1

          do layer=1,num_layers
            ! compartment 2 - process 1: incomplete division
            prop2(1, layer, x1, x2) = up_to_half(lt2(layer, x1, x2)) * div_2
            ! compartment 2 - process 2: fragmentation
            prop2(2, layer, x1, x2) = larger_one(lt2(layer, x1, x2)) * (lt2(layer, x1, x2) - 1) * frag_2
            ! compartment 2 - process 3: dedifferentiation
            prop2(3, layer, x1, x2) = larger_zero(lt2(layer, x1, x2)) * dediff_2
            ! compartment 2 - process 4: differentiation
            prop2(4, layer, x1, x2) = larger_zero(lt2(layer, x1, x2)) * diff_2

            ! compartment 3 - process 1: incomplete division
            prop3(1, layer, x1, x2) = up_to_half(lt3(layer, x1, x2)) * div_3
            ! compartment 3 - process 2: fragmentation - proportional to the number of bridges
            prop3(2, layer, x1, x2) = larger_one(lt3(layer, x1, x2)) * (lt3(layer, x1, x2) - 1) * frag_3
            ! compartment 3 - process 3: differentiation
            prop3(3, layer, x1, x2) = larger_zero(lt3(layer, x1, x2)) * diff_3
            ! compartment 3 - process 4: dedifferentiation
            prop3(4, layer, x1, x2) = larger_zero(lt3(layer, x1, x2)) * dediff_3
          end do
        end do
      end do

      ! -- advance time and find reaction
      total_prop = prop_osc + sum(prop1) + sum(prop2) + sum(prop3)

      call random_number(time_step)
      time_step = - (1. / total_prop) * log(time_step)
      t = t + time_step

      call random_number(reaction)
      reaction = total_prop * reaction

      ! -- reactions
      p_sum = 0.
      searching = .true.

      x1 = 0
      do while(x1 < L1.and.searching)
        x1 = x1 + 1
        x2 = 0
        do while(x2 < L2.and.searching)
          x2 = x2 + 1

          ! === COMPARTMENT 1 =====

          ! compartment 1 - process 2: fragmentation
          p_sum = p_sum + prop1(2, x1, x2)
          if (reaction < p_sum.and.searching) then
            ! remember the original state of the lattice site
            last_bp = 0
            orig_length = lt1(x1, x2)

            ! loop through all bridges
            do i = 1, orig_length - 1
              ! determine whether bridge breaks
              call random_number(die)
              if (die < 0.5) then
                ! truncate original unit
                lt1(x1, x2) = orig_length - i

                ! choose neighboring lt site for invasion
                call random_neighbor(L_mig, x1, x2, L1, L2, y1, y2)

                ! search for free layer in upper lattice
                inv = free_layer(lt2(:, y1, y2), num_layers)

                ! differentiate invaded site
                lt2(inv, y1, y2) = lt1(y1, y2)
                pm2(inv, y1, y2) = pm1(y1, y2)
                id2(inv, y1, y2) = id1(y1, y2)

                ! track flux: X -> Y
                if (equilibrated) then
                  cell_flux(1, 2) = cell_flux(1, 2) + lt1(y1, y2)
                end if

                ! invade neighboring site
                lt1(y1, y2) = i - last_bp
                pm1(y1, y2) = pm1(x1, x2)
                id1(y1, y2) = id1(x1, x2)

                ! track last breakpoint
                last_bp = i
              end if
            end do

            searching = .false.
            exit
          end if

          do layer = 1, num_layers
            ! === COMPARTMENT 2 ===

            ! compartment 2 - process 1: incomplete division
            p_sum = p_sum + prop2(1, layer, x1, x2)
            if (reaction < p_sum.and.searching) then
              ! track added number of cells: Y -> 2Y
              if (equilibrated) then
                cell_flux(2, 2) = cell_flux(2, 2) + lt2(layer, x1, x2)
              end if

              total_div_events = total_div_events + lt2(layer, x1, x2)
              lt2(layer, x1, x2) = lt2(layer, x1, x2) * 2
              pm2(layer, x1, x2) = pm2(layer, x1, x2) + 1

              searching = .false.
              exit
            end if

            ! compartment 2 - process 2: fragmentation
            p_sum = p_sum + prop2(2, layer, x1, x2)
            if (reaction < p_sum.and.searching) then
              ! remember the original state of the lattice site
              last_bp = 0
              orig_length = lt2(layer, x1, x2)

              ! loop through all bridges
              do i = 1, orig_length - 1
                ! determine whether bridge breaks
                call random_number(die)
                if (die < 0.5) then
                  ! truncate original unit
                  lt2(layer, x1, x2) = orig_length - i

                  ! choose neighboring lt site for invasion
                  call random_neighbor(L_mig, x1, x2, L1, L2, y1, y2)

                  ! search for free layer
                  inv = free_layer(lt2(:, y1, y2), num_layers)

                  ! invade neighboring site
                  lt2(inv, y1, y2) = i - last_bp
                  pm2(inv, y1, y2) = pm2(layer, x1, x2)
                  id2(inv, y1, y2) = id2(layer, x1, x2)

                  ! track last breakpoint
                  last_bp = i
                end if
              end do

              searching = .false.
              exit
            end if

            ! compartment 2 - process 3: dedifferentiation
            p_sum = p_sum + prop2(3, layer, x1, x2)
            if (reaction < p_sum.and.searching) then
              ! differentiate site on lt1: swap with dedifferentiating unit
              temp_lt = lt2(layer, x1, x2)
              temp_pm = pm2(layer, x1, x2)
              temp_id = id2(layer, x1, x2)

              lt2(layer, x1, x2) = lt1(x1, x2)
              pm2(layer, x1, x2) = pm1(x1, x2)
              id2(layer, x1, x2) = id1(x1, x2)

              ! track flux: X -> Y
              if (equilibrated) then
                cell_flux(1, 2) = cell_flux(1, 2) + lt1(x1, x2)
              end if

              ! dedifferentiate
              lt1(x1, x2) = temp_lt
              pm1(x1, x2) = temp_pm
              id1(x1, x2) = temp_id

              ! track flux: Y -> X
              if (equilibrated) then
                cell_flux(2, 1) = cell_flux(2, 1) + temp_lt
              end if

              searching = .false.
              exit
            end if

            ! compartment 2 - process 4: differentiation
            p_sum = p_sum + prop2(4, layer, x1, x2)
            if (reaction < p_sum.and.searching) then
              ! search for free layer
              inv = free_layer(lt3(:, x1, x2), num_layers)

              ! differentiate site on lt2
              lt3(inv, x1, x2) = lt2(layer, x1, x2)
              pm3(inv, x1, x2) = pm2(layer, x1, x2)
              id3(inv, x1, x2) = id2(layer, x1, x2)

              ! track flux: Y -> Z
              if (equilibrated) then
                cell_flux(2, 3) = cell_flux(2, 3) + lt2(layer, x1, x2)
              end if

              lt2(layer, x1, x2) = 0
              pm2(layer, x1, x2) = 0
              id2(layer, x1, x2) = 0

              searching = .false.
              exit
            end if

            ! === COMPARTMENT 3 ===

            ! compartment 3 - process 1: incomplete division
            p_sum = p_sum + prop3(1, layer, x1, x2)
            if (reaction < p_sum.and.searching) then
              ! track added number of cells: Z -> 2Z
              if (equilibrated) then
                cell_flux(3, 3) = cell_flux(3, 3) + lt3(layer, x1, x2)
              end if

              total_div_events = total_div_events + lt3(layer, x1, x2)
              lt3(layer, x1, x2) = lt3(layer, x1, x2) * 2
              pm3(layer, x1, x2) = pm3(layer, x1, x2) + 1

              searching = .false.
              exit
            end if

            ! compartment 3 - process 2: fragmentation
            p_sum = p_sum + prop3(2, layer, x1, x2)
            if (reaction < p_sum.and.searching) then
              ! remember the original state of the lattice site
              last_bp = 0
              orig_length = lt3(layer, x1, x2)

              ! loop through all bridges
              do i = 1, orig_length - 1
                ! determine whether bridge breaks
                call random_number(die)
                if (die < 0.5) then
                  ! truncate original unit
                  lt3(layer, x1, x2) = orig_length - i

                  ! choose neighboring lt site for invasion
                  call random_neighbor(L_mig, x1, x2, L1, L2, y1, y2)

                  ! search for free layer
                  inv = free_layer(lt3(:, y1, y2), num_layers)

                  ! invade neighboring site
                  lt3(inv, y1, y2) = i - last_bp
                  pm3(inv, y1, y2) = pm3(layer, x1, x2)
                  id3(inv, y1, y2) = id3(layer, x1, x2)

                  ! track last breakpoint
                  last_bp = i
                end if
              end do

              searching = .false.
              exit
            end if

            ! compartment 3 - process 3: differentiation
            p_sum = p_sum + prop3(3, layer, x1, x2)
            if (reaction < p_sum.and.searching) then
              ! track added number of cells: Z -> (/)
              if (equilibrated) then
                cell_flux(3, 4) = cell_flux(3, 4) + lt3(layer, x1, x2)
              end if

              ! track number of divisions upon exit
              total_n_div_exit = total_n_div_exit + pm3(layer, x1, x2) * lt3(layer, x1, x2)
              lost_cells = lost_cells + lt3(layer, x1, x2)

              lt3(layer, x1, x2) = 0
              pm3(layer, x1, x2) = 0
              id3(layer, x1, x2) = 0

              searching = .false.
              exit
            end if

            ! compartment 3 - process 4: dedifferentiation
            p_sum = p_sum + prop3(4, layer, x1, x2)
            if (reaction < p_sum.and.searching) then
              ! track added number of cells: Z -> (/)
              if (equilibrated) then
                cell_flux(3, 2) = cell_flux(3, 2) + lt3(layer, x1, x2)
              end if

              ! search for free layer
              inv = free_layer(lt2(:, x1, x2), num_layers)

              lt2(inv, x1, x2) = lt3(layer, x1, x2)
              pm2(inv, x1, x2) = pm3(layer, x1, x2)
              id2(inv, x1, x2) = id3(layer, x1, x2)

              lt3(layer, x1, x2) = 0
              pm3(layer, x1, x2) = 0
              id3(layer, x1, x2) = 0

              searching = .false.
              exit
            end if
          end do
        end do
        if (.not.searching) exit
      end do

      ! advance oscillator
      p_sum = p_sum + prop_osc
      if (reaction < p_sum.and.searching) then
        osc = osc + 1

        x1 = 0
        do while(x1<L1.and.searching)
          x1 = x1 + 1
          x2 = 0
          do while(x2<L2.and.searching)
            x2 = x2 + 1

            if ((lt1(x1, x2) <= num_syncytia_2).and.(modulo(osc - nint(real(m,8) * real(x2,8) / real(L2,8)), m) == 0)) then
              total_div_events = total_div_events + lt1(x1, x2)
              ! determine whether unit differentiates
              call random_number(die)
              if ((die < prob_diff_1).and.(sum(lt2(:,x1,x2)) > 0)) then
                ! find occupied layer
                inv = occ_layer(lt2(:, x1, x2), num_layers)

                ! track added number of cells: X -> 2Y
                ! (counted as X -> 2X followed by 2X -> 2Y)
                if (equilibrated) then
                  cell_flux(1, 1) = cell_flux(1, 1) + lt1(x1, x2)
                  cell_flux(1, 2) = cell_flux(1, 2) + lt1(x1, x2) * 2
                  cell_flux(2, 1) = cell_flux(2, 1) + lt2(inv, x1, x2)
                end if

                temp_lt = lt2(inv, x1, x2)
                temp_pm = pm2(inv, x1, x2)
                temp_id = id2(inv, x1, x2)

                lt2(inv, x1, x2) = lt1(x1, x2) * 2
                pm2(inv, x1, x2) = pm1(x1, x2) + 1
                id2(inv, x1, x2) = id1(x1, x2)

                lt1(x1, x2) = temp_lt
                pm1(x1, x2) = temp_pm
                id1(x1, x2) = temp_id

                ! count cell fate as "wanted to differentiate and did"
                if (equilibrated) then
                  total_dupl_fates(2) = total_dupl_fates(2) + 1
                end if
              else
                ! track added number of cells: X -> 2X
                if (equilibrated) then
                  cell_flux(1,1) = cell_flux(1,1) + lt1(x1, x2)
                end if

                lt1(x1, x2) = lt1(x1, x2) * 2
                pm1(x1, x2) = pm1(x1, x2) + 1

                if (equilibrated) then
                  if (die < prob_diff_1) then
                    ! count cell fate as "wanted to differentiate but couldn't"
                    total_dupl_fates(3) = total_dupl_fates(3) + 1
                  else
                    ! count cell fate as "didn't want to differentiate"
                    total_dupl_fates(1) = total_dupl_fates(1) + 1
                  end if
                end if
              end if
            end if
          end do
        end do

        searching = .false.
      end if

      ! --- data storage ---

      rp = floor(t / dt)
      if (rp > lrp.and.rp < num_rps.and.equilibrated) then
        ! - add to global distributions
        do rpi = lrp + 1, rp
          avg_csd(rpi, :, :, :, :) = avg_csd(rpi, :, :, :, :) + real(last_csd(:, :, :, :), 8)
          avg_cmp(rpi, :, :) = avg_cmp(rpi, :, :) + last_cmp(:, :)
          avg_labelled_cmp(rpi, :, :) = avg_labelled_cmp(rpi, :, :) + last_num_labelled_units(:, :)
        end do

        avg_lf2(lrp+1:rp) = avg_lf2(lrp+1:rp) + last_lf2 / real(num_sites, 8)
        avg_lf3(lrp+1:rp) = avg_lf3(lrp+1:rp) + last_lf3 / real(num_sites, 8)

        do rpi = lrp + 1, rp
          avg_lu(rpi, :) = avg_lu(rpi, :) + last_lu(:)
          avg_lc_frac(rpi, :) = avg_lc_frac(rpi, :) + last_lc_frac(:)
        end do

        do c = 1, 3
          do s = 1, num_syncytia
            if (last_num_units(c,s)>0) then
              nrm_marked_ratio(lrp+1:rp, c, s) = nrm_marked_ratio(lrp+1:rp, c, s) + 1
              avg_marked_ratio(lrp+1:rp, c, s) = avg_marked_ratio(lrp+1:rp, c, s) + real(last_num_marked_units(c, s), 8) / real(last_num_units(c, s), 8)
            end if
          end do
        end do

        avg_lost_ids(lrp+1:rp) = avg_lost_ids(lrp+1:rp) + last_avg_lost_ids

        do rpi = lrp + 1, rp
          avg_osc(rpi) = avg_osc(rpi) + last_osc
          avg_n_div(rpi, :) = avg_n_div(rpi, :) + last_n_div(:)
          avg_lost_cells(rpi) = avg_lost_cells(rpi) + last_lost_cells
          avg_total_div_events(rpi) = avg_total_div_events(rpi) + last_total_div_events
          avg_total_n_div_exit(rpi) = avg_total_n_div_exit(rpi) + last_total_n_div_exit
          avg_n_div_dist(rpi, :, :) = avg_n_div_dist(rpi, :, :) + last_n_div_dist(:, :)
#ifdef labelcomp
          avg_origin(rpi, :) = avg_origin(rpi, :) + last_origin(:)
#endif
          avg_cell_flux(rpi, :, :) = avg_cell_flux(rpi, :, :) + cell_flux(:, :)
          avg_div_total(rpi, :) = avg_div_total(rpi, :) + last_div_total(:)
        end do

        lrp = rp
      end if

      ! --- data analysis ---

      ! obtain average filling of the second c and ratio of divided cells
      lf2 = 0
      lf3 = 0
      num_units = 0
      num_marked_units = 0
      num_labelled_units = 0
      lu = 0
      lc = 0
      last_n_div_dist = 0
      origin = 0
      do x1 = 1, L1
        do x2 = 1, L2
          ! compartment 1
          c = 1
          unit_type = lt1(x1, x2)
          if (unit_type > 0) then
            num_units(c, unit_type) = num_units(c, unit_type) + 1
            if (pm1(x1, x2) > 0) then
              num_marked_units(c, unit_type) = num_marked_units(c, unit_type) + 1
            end if
            if (id1(x1, x2)>0) then
              lu(c) = lu(c) + 1
              lc(c) = lc(c) + unit_type
              num_labelled_units(c, unit_type) = num_labelled_units(c, unit_type) + 1
            end if
            last_n_div_dist(c, pm1(x1, x2)) = last_n_div_dist(c, pm1(x1, x2)) + unit_type
          end if


#ifdef labelcomp
          ! origin of cells: only of labelling mode is 'labelcomp'
          origin(id1(x1, x2)) = origin(id1(x1, x2)) + lt1(x1, x2)
#endif

          do layer = 1, num_layers
                ! compartment 2
            c = 2
            unit_type = lt2(layer, x1, x2)
            lf2 = lf2 + real(unitstep(unit_type), 8)
            if (unit_type > 0) then
              num_units(c, unit_type) = num_units(c, unit_type) + 1
              if (pm2(layer, x1, x2) > 0) then
                num_marked_units(c,unit_type) = num_marked_units(c,unit_type) + 1
              end if
              if (id2(layer, x1, x2) > 0) then
                lu(c) = lu(c) + 1
                lc(c) = lc(c) + unit_type
                num_labelled_units(c,unit_type) = num_labelled_units(c,unit_type) + 1
              end if
              last_n_div_dist(c, pm2(layer, x1, x2)) = last_n_div_dist(c, pm2(layer, x1, x2)) + unit_type
            end if

            ! compartment 3
            c = 3
            unit_type = lt3(layer, x1, x2)
            lf3 = lf3 + real(unitstep(unit_type), 8)
            if (unit_type > 0) then
              num_units(c, unit_type) = num_units(c, unit_type) + 1
              if (pm3(layer, x1, x2) > 0) then
                num_marked_units(c, unit_type) = num_marked_units(c, unit_type) + 1
              end if
              if (id3(layer, x1, x2) > 0) then
                lu(c) = lu(c) + 1
                lc(c) = lc(c) + unit_type
                num_labelled_units(c, unit_type) = num_labelled_units(c, unit_type) + 1
              end if
              last_n_div_dist(c, pm3(layer, x1, x2)) = last_n_div_dist(c, pm3(layer, x1, x2)) + unit_type
            end if

            ! origin of cells: only of labelling mode is 'labelcomp'
#ifdef labelcomp
            origin(id2(layer, x1, x2)) = origin(id2(layer, x1, x2)) + lt2(layer, x1, x2)
            origin(id3(layer, x1, x2)) = origin(id3(layer, x1, x2)) + lt3(layer, x1, x2)
#endif
          end do
        end do
      end do

      ! - determine clone sizes and cmps for each barcode
      call analyze_clones(max_cs, num_ids, num_layers, num_syncytia, L1, L2, lt1, lt2, lt3, id1, id2, id3, csd, cmp)

      last_osc = osc
      last_csd = csd
      last_avg_lost_ids = real(csd(0, 0, 0, UNITS),8) / real(num_ids, 8)
      last_cmp = cmp
      last_lf2 = lf2
      last_lf3 = lf3
      last_lu = lu
      last_div_total(1) = real(sum(pm1), 8)
      last_div_total(2) = real(sum(pm2), 8)
      last_div_total(3) = real(sum(pm3), 8)
      if (sum(lt1) > 0) then
        last_lc_frac(1) = real(lc(1), 8) / real(sum(lt1), 8)
        last_n_div(1) = real(sum(pm1 * lt1), 8) / real(sum(lt1), 8)
      else
        last_lc_frac(1) = real(lc(1), 8)
        last_n_div(1) = real(sum(pm1 * lt1), 8)
      end if
      if (sum(lt2) > 0) then
        last_lc_frac(2) = real(lc(2), 8) / real(sum(lt2), 8)
        last_n_div(2) = real(sum(pm2 * lt2), 8) / real(sum(lt2), 8)
      else
        last_lc_frac(2) = real(lc(2), 8)
        last_n_div(2) = real(sum(pm2 * lt2), 8)
      end if
      if (sum(lt3) > 0) then
        last_lc_frac(3) = real(lc(3), 8) / real(sum(lt3), 8)
        last_n_div(3) = real(sum(pm3 * lt3), 8) / real(sum(lt3), 8)
      else
        last_lc_frac(3) = real(lc(3), 8)
        last_n_div(3) = real(sum(pm3 * lt3), 8)
      end if
      last_num_marked_units = num_marked_units
      last_num_labelled_units = num_labelled_units
      last_num_units = num_units
      last_lost_cells = lost_cells
      last_total_div_events = total_div_events
      last_total_n_div_exit = total_n_div_exit
      last_origin = origin
      last_cell_flux = cell_flux
    end do

  ! end realizations
  end do

  ! --- data analysis ---

  ! normalise global quantities
  avg_osc = avg_osc / real(nr, 8)
  avg_csd = avg_csd / real(nr, 8)
  avg_cmp = avg_cmp / real(nr, 8)
  avg_lf2 = avg_lf2 / real(nr, 8)
  avg_lf3 = avg_lf3 / real(nr, 8)
  avg_lu = avg_lu / real(nr, 8)
  avg_lc_frac = avg_lc_frac / real(nr, 8)
  avg_lost_ids = avg_lost_ids / real(nr, 8)
  avg_n_div = avg_n_div / real(nr, 8)
  avg_div_total = avg_div_total / real(nr, 8)
  avg_lost_cells = avg_lost_cells / (real(nr, 8) * num_sites)
  avg_total_div_events = avg_total_div_events / (real(nr, 8) * num_sites)
  avg_total_n_div_exit = avg_total_n_div_exit / (real(nr, 8) * num_sites)
  avg_n_div_dist = avg_n_div_dist / real(nr, 8)
  avg_origin = avg_origin / real(nr, 8)
  avg_cell_flux = avg_cell_flux / real(nr, 8)
  avg_labelled_cmp = avg_labelled_cmp / real(nr, 8)

  do rp=1,num_rps
  do c=1,3
    do s=1,num_syncytia
    if (nrm_marked_ratio(rp, c, s) > 0) then
      avg_marked_ratio(rp, c, s) = avg_marked_ratio(rp, c, s) / real(nrm_marked_ratio(rp, c, s), 8)
    else
      avg_marked_ratio(rp, c, s) = -1.
    end if
    end do
  end do
  end do

  ! write parameter file
  open(unit=1, file=trim(directory) // 'parameters.csv')
    call write_parameters(1, num_layers, num_syncytia, nu, m, frag_1, div_2, frag_2, prob_dediff, conv_2, div_3, frag_3, diff_3, dediff_3, L1, L2, dt, sim_time, eq_time, nr, seed)
  close(1)

  ! --- data storage ---

#ifdef single
  ! write average clone size distribution
  do k = 0, num_rps
    ! --- X (rows) vs. Y + Z (columns)
    write(filename, '(a,i5.5,a)') trim(directory) // 'avg_clone_sizes_units_plvap_', k, '.csv'
    open(unit=1,file=trim(filename))
    write(filename, '(a,i5.5,a)') trim(directory) // 'avg_clone_sizes_cells_plvap_', k, '.csv'
    open(unit=2,file=trim(filename))

    do n = 0, max_cs
      do unit_type = 1, 2
        ! marginalise to obtain an X+/X- distribution
        marg_avg_csd = 0.
        do n1 = 0, max_cs
          do n2 = 0, max_cs
            nt = n1 + n2
            if (nt <= max_cs) then
              marg_avg_csd(nt) = marg_avg_csd(nt) + avg_csd(k, n, n1, n2, unit_type)
            end if
          end do
        end do

        write(unit_type, *) marg_avg_csd
      end do
    end do

    ! --- Z (rows) vs. X + Y (columns)
    write(filename, '(a,i5.5,a)') trim(directory) // 'avg_clone_sizes_units_sox3_', k, '.csv'
    open(unit=1,file=trim(filename))
    write(filename, '(a,i5.5,a)') trim(directory) // 'avg_clone_sizes_cells_sox3_', k, '.csv'
    open(unit=2,file=trim(filename))

    do n = 0, max_cs
      do unit_type = 1, 2
        ! marginalise to obtain an X+/X- distribution
        marg_avg_csd = 0.
        do n1 = 0, max_cs
          do n2 = 0, max_cs
            nt = n1 + n2
            if (nt <= max_cs) then
              marg_avg_csd(nt) = marg_avg_csd(nt) + avg_csd(k, n1, n2, n, unit_type)
            end if
          end do
        end do

        write(unit_type, *) marg_avg_csd
      end do
    end do

    close(1)
    close(2)
  end do
#endif

  ! write average cmp and ratios of divided cells
  open(unit=1,file=trim(directory) // 'avg_composition.csv')
  open(unit=2,file=trim(directory) // 'avg_marked_ratio.csv')
  open(unit=3,file=trim(directory) // 'avg_labelled_composition.csv')
    do k = 0, num_rps
      write(3, '(f17.6,a1)', advance='no') k * dt, ','
      do c = 1, 3
        do s = 1, num_syncytia
          if ((c == 3).and.(s == num_syncytia)) then
            write(1, '(f17.6)') avg_cmp(k, c, s)
            write(2, '(f17.6)') avg_marked_ratio(k, c, s)
            write(3, '(f17.6)') avg_labelled_cmp(k, c, s)
          else
            write(1, '(f17.6,a1)', advance='no') avg_cmp(k, c, s), ','
            write(2, '(f17.6,a1)', advance='no') avg_marked_ratio(k, c, s), ','
            write(3, '(f17.6,a1)', advance='no') avg_labelled_cmp(k, c, s), ','
          end if
        end do
      end do
    end do
  close(1)
  close(2)
  close(3)

#ifndef labelcomp
  do uc = 1, 2
    if (uc == 1) then
      open(unit=1, file=trim(directory) // 'avg_clone_content_total_units.csv')
      open(unit=2, file=trim(directory) // 'avg_clone_content_surv_units.csv')
    end if
    if (uc == 2) then
      open(unit=1, file=trim(directory) // 'avg_clone_content_total_cells.csv')
      open(unit=2, file=trim(directory) // 'avg_clone_content_surv_cells.csv')
    end if

    do k = 0, num_rps
      ! -- total population
      curdist = avg_csd(k,:,:,:,uc)
      if (sum(curdist) > 0.) curdist = curdist / sum(curdist)
      ! compute average clone content
      avg_content = 0.
      do n = 0, max_cs
        avg_content(1) = avg_content(1) + real(n, 8) * sum(curdist(n, :, :))
        avg_content(2) = avg_content(2) + real(n, 8) * sum(curdist(:, n, :))
        avg_content(3) = avg_content(3) + real(n, 8) * sum(curdist(:, :, n))
      end do
      write(1, '(f17.6,a1,f17.6,a1,f17.6,a1,f17.6)') k * dt, ',', avg_content(1), ',', avg_content(2), ',', avg_content(3)

      ! -- surviving population
      ! normalise clone size distribution to surviving population
      curdist = avg_csd(k, :, :, :, uc)
      curdist(0, 0, 0) = 0.
      if (sum(curdist) > 0.) curdist = curdist / sum(curdist)

      ! compute average clone content
      avg_content = 0.
      do n = 0, max_cs
        avg_content(1) = avg_content(1) + real(n, 8) * sum(curdist(n, :, :))
        avg_content(2) = avg_content(2) + real(n, 8) * sum(curdist(:, n, :))
        avg_content(3) = avg_content(3) + real(n, 8) * sum(curdist(:, :, n))
      end do
      write(2,'(f17.6,a1,f17.6,a1,f17.6,a1,f17.6)') k * dt, ',', avg_content(1), ',', avg_content(2), ',', avg_content(3)
    end do
  end do
#endif

  open(unit=1, file=trim(directory) // 'avg_osc.csv')
    do k = 0, num_rps
      write(1, '(f17.6, f17.6)') k * dt, avg_osc(k)
    end do
  close(1)

  open(unit=1, file=trim(directory) // 'avg_labelled_units.csv')
    do k = 0, num_rps
      write(1, '(f17.6, f17.6, f17.6, f17.6)') k * dt, avg_lu(k, 1), avg_lu(k, 2), avg_lu(k, 3)
    end do
  close(1)

  open(unit=1, file=trim(directory) // 'avg_labelled_cells_frac.csv')
    do k = 0, num_rps
      write(1, '(f17.6, f17.6, f17.6, f17.6)') k * dt, avg_lc_frac(k, 1), avg_lc_frac(k, 2), avg_lc_frac(k, 3)
    end do
  close(1)

#ifndef labelcomp
  open(unit=1,file=trim(directory) // 'avg_lost_ids.csv')
    do k = 0, num_rps
      write(1, '(f17.6, f17.6)') k * dt, avg_lost_ids(k)
    end do
  close(1)
#endif

  open(unit=1, file=trim(directory) // 'avg_diff_characteristics.csv')
    do k = 0, num_rps
      write(1, '(f17.6, f17.6, f17.6, f17.6)') k * dt, avg_lost_cells(k), avg_total_n_div_exit(k), avg_total_div_events(k)
    end do
  close(1)

  open(unit=1,file=trim(directory) // 'avg_n_div.csv')
    do k = 0, num_rps
      write(1, '(f17.6, f17.6, f17.6, f17.6)') k * dt, avg_n_div(k, 1), avg_n_div(k, 2), avg_n_div(k, 3)
    end do
  close(1)

  open(unit=1,file=trim(directory) // 'avg_div_total.csv')
    do k = 0, num_rps
      write(1, '(f17.6, f17.6, f17.6, f17.6)') k * dt, avg_div_total(k, 1), avg_div_total(k, 2), avg_div_total(k, 3)
    end do
  close(1)

#ifdef labelcomp
  open(unit=1, file=trim(directory) // 'avg_origin.csv')
    do k = 0, num_rps
      write(1, '(f17.6, f17.6, f17.6, f17.6)') k * dt, avg_origin(k, 1), avg_origin(k, 2), avg_origin(k, 3)
    end do
  close(1)
#endif

  open(unit=1,file=trim(directory) // 'avg_cell_flux.csv')
    do k = 0, num_rps
      write(1, '(f17.6, f17.6, f17.6, f17.6, f17.6, f17.6, f17.6, f17.6, f17.6)') k * dt, avg_cell_flux(k, 1, 1), avg_cell_flux(k, 1, 2), avg_cell_flux(k, 2, 1), avg_cell_flux(k, 2, 2), avg_cell_flux(k, 2, 3), avg_cell_flux(k, 3, 2), avg_cell_flux(k, 3, 3), avg_cell_flux(k, 3, 4)
    end do
  close(1)

  c = 1
  open(unit=1,file=trim(directory) // 'avg_n_div_dist_cells_x.csv')
    do k = 0, num_rps
      write(1, '(f17.6,a1)', advance='no') k * dt, '  '
      do i = 0, max_n_div
        write(1, '(f17.6,a1)', advance='no') avg_n_div_dist(k, c, i), '  '
      end do
      write(1, *)
    end do
  close(1)

  c = 2
  open(unit=1,file=trim(directory) // 'avg_n_div_dist_cells_y.csv')
    do k = 0, num_rps
      write(1, '(f17.6,a1)', advance='no') k * dt, '  '
      do i = 0, max_n_div
        write(1, '(f17.6,a1)', advance='no') avg_n_div_dist(k, c, i), '  '
      end do
      write(1, *)
    end do
  close(1)

  c = 3
  open(unit=1,file=trim(directory) // 'avg_n_div_dist_cells_z.csv')
    do k = 0, num_rps
      write(1, '(f17.6,a1)', advance='no') k * dt, '  '
      do i = 0, max_n_div
        write(1, '(f17.6,a1)', advance='no') avg_n_div_dist(k, c, i), '  '
      end do
      write(1, *)
    end do
  close(1)

  ! write global quantities
  open(unit=1, file=trim(directory) // 'avg_layer_filling.csv')
    write(1, *) avg_lf2
    write(1, *) avg_lf3
  close(1)

  open(unit=1, file=trim(directory) // 'total_dupl_fates.csv')
    write(1,'(i10,a1,i10,a1,i10)') total_dupl_fates(1), ',', total_dupl_fates(2), ',', total_dupl_fates(3)
  close(1)

end program germline


subroutine write_parameters(output_unit, num_layers, num_syncytia, nu, m, frag_1, div_2, frag_2, prob_dediff, conv_2, div_3, frag_3, diff_3, dediff_3, L1, L2, dt, sim_time, eq_time, nr, seed)
  implicit none

  integer, intent(in) :: output_unit
  real(kind=8), intent(in) :: nu, frag_1, div_2, frag_2, prob_dediff, conv_2, div_3, frag_3, diff_3, dediff_3, dt, sim_time, eq_time
  integer, intent(in) :: num_layers, num_syncytia, L1, L2, nr, m
  integer, intent(in) :: seed

  write(output_unit,'(a20,f17.6)') 'nu, ', nu
  write(output_unit,'(a20,i15)')   'm, ', m
  write(output_unit,'(a20,f17.6)') 'frag_1, ', frag_1
  write(output_unit,'(a20,f17.6)') 'div_2, ', div_2
  write(output_unit,'(a20,f17.6)') 'frag_2, ', frag_2
  write(output_unit,'(a20,f17.6)') 'prob_dediff, ', prob_dediff
  write(output_unit,'(a20,f17.6)') 'conv_2, ', conv_2
  write(output_unit,'(a20,f17.6)') 'div_3, ', div_3
  write(output_unit,'(a20,f17.6)') 'frag_3, ', frag_3
  write(output_unit,'(a20,f17.6)') 'diff_3, ', diff_3
  write(output_unit,'(a20,f17.6)') 'dediff_3, ', dediff_3

  write(output_unit,'(a20,i15)')   'L1, ', L1
  write(output_unit,'(a20,i15)')   'L2, ', L2

  write(output_unit,'(a20,f17.6)') 'delta_t, ', dt
  write(output_unit,'(a20,f17.6)') 'simulation_time, ', sim_time
  write(output_unit,'(a20,f17.6)') 'equilibration_time, ', eq_time
  write(output_unit,'(a20,i15)')   'num_realizations, ', nr
  write(output_unit,'(a20,i15)')   'seed, ', seed

  write(output_unit,'(a20,i15)')   'num_layers, ', num_layers
  write(output_unit,'(a20,i15)')   'num_syncytia, ', num_syncytia
end subroutine write_parameters


subroutine analyze_clones(max_cs, num_ids, num_layers, num_syncytia, L1, L2, lt1, lt2, lt3, id1, id2, id3, csd, cmp)
  implicit none

  integer, parameter :: UNITS = 1, CELLS = 2

  integer, intent(in) :: max_cs, num_ids, num_layers, num_syncytia, L1, L2
  integer, dimension(1:L1, 1:L2), intent(in) :: lt1, id1
  integer, dimension(1:num_layers, 1:L1, 1:L2), intent(in) :: lt2, lt3, id2, id3

  integer :: x1, x2, layer, c, a, unit_type
  integer, dimension(0:num_ids, 1:2, 1:3) :: clone_size
  integer, dimension(1:3) :: cs

  integer, dimension(0:max_cs, 0:max_cs, 0:max_cs, 1:2), intent(out) :: csd
  integer, dimension(1:3, 0:num_syncytia), intent(out) :: cmp

  integer :: unitstep

  ! determine sizes of individual clones
  clone_size = 0
  cmp = 0

  do x1 = 1, L1
    do x2 = 1, L2
      a = id1(x1, x2)

      clone_size(a, UNITS, 1) = clone_size(a, UNITS, 1) + unitstep(lt1(x1, x2))
      clone_size(a, CELLS, 1) = clone_size(a, CELLS, 1) + lt1(x1, x2)
      cmp(1, lt1(x1, x2)) = cmp(1, lt1(x1, x2)) + 1

      do layer = 1, num_layers
        c = 2
        a = id2(layer, x1, x2)
        clone_size(a, UNITS, c) = clone_size(a, UNITS, c) + unitstep(lt2(layer, x1, x2))
        clone_size(a, CELLS, c) = clone_size(a, CELLS, c) + lt2(layer, x1, x2)
        cmp(c, lt2(layer, x1, x2)) = cmp(c, lt2(layer, x1, x2)) + 1

        c = 3
        a = id3(layer, x1, x2)
        clone_size(a, UNITS, c) = clone_size(a, UNITS, c) + unitstep(lt3(layer, x1, x2))
        clone_size(a, CELLS, c) = clone_size(a, CELLS, c) + lt3(layer, x1, x2)
        cmp(c, lt3(layer, x1, x2)) = cmp(c, lt3(layer, x1, x2)) + 1
      end do
    end do
  end do

  ! determine clone size distribution
  csd = 0

  do a = 1, num_ids
    do unit_type = 1, 2
      cs(:) = clone_size(a, unit_type, :)
      if ((cs(1) <= max_cs).and.(cs(2) <= max_cs).and.(cs(3) <= max_cs)) then
        csd(cs(1), cs(2), cs(3), unit_type) = csd(cs(1), cs(2), cs(3), unit_type) + 1
      end if
    end do
  end do
end subroutine


subroutine write_to_files(num_layers, num_sites, num_syncytia, L1, L2, lt, id, clone_size_dist, cmp, rp, directory_r)
  implicit none

  integer, intent(in) :: num_layers, num_sites, num_syncytia, L1, L2, rp
  integer, dimension(1:num_layers, 1:L1, 1:L2), intent(in)  :: lt, id
  integer,  dimension(1:num_layers, 0:num_sites*num_layers), intent(in) :: clone_size_dist
  integer, dimension(1:num_layers, 0:num_syncytia), intent(in) :: cmp
  character(255), intent(in) :: directory_r

  integer :: x1, x2, layer, s, n
  character(255) :: filename_frame

  ! - cmp
  do layer = 1, num_layers
    do s = 1, num_syncytia
      write(2, '(i10)', advance='no') cmp(layer, s)
    end do
  end do
  write(2, *)

  ! - clone size distributions
  do layer = 1, num_layers
    do n = 0, num_sites
      write(2+layer, '(i10)', advance='no') clone_size_dist(layer, n)
    end do
    write(2+layer, *)
  end do

  ! write entire lt configuration to files
  write(filename_frame, '(a,i5.5,a)') trim(directory_r) // 'frame_', rp, '.csv'
  open(unit=1, file=trim(filename_frame))
    ! write cells
    do layer = 1, num_layers
      do x2 = 1, L2
        write(1, *) lt(layer, :, x2)
      end do
      write(1, *)
    end do

    ! write clone IDs
    do layer = 1, num_layers
      do x2 = 1, L2
        write(1, *) id(layer, :, x2)
      end do
      write(1, *)
    end do
  close(1)
end subroutine


! note: modulo which starts from 1 instead of 0 is given by rmod(n,m) = mod(n-1,m)+1
function rmod(n, m)
  implicit none

  integer :: rmod
  integer, intent(in) :: n, m

  rmod = modulo(n-1,m)+1
  return
end function


function unitstep(x)
  implicit none

  integer :: unitstep
  integer :: x

  unitstep = 0
  if (x > 0) unitstep = 1
  return
end function


function free_layer(array, length)
  implicit none

  integer :: free_layer

  integer, intent(in) :: length
  integer, dimension(1:length), intent(in) :: array

  integer :: i
  real(kind=8) :: die

  free_layer = -1
  do i = 1, length
    if (array(i) == 0) then
      free_layer = i
      exit
    end if
  end do
  if (free_layer == -1) then
    call random_number(die)
    free_layer = 1 + nint(real(length - 1, 8) * die)
  end if
  return
end function


function occ_layer(array, length)
  implicit none

  integer :: occ_layer

  integer, intent(in) :: length
  integer, dimension(1:length), intent(in) :: array

  integer, dimension(1:length) :: candidates

  integer :: i, c
  real(kind=8) :: die

  c = 0
  do i = 1, length
    if (array(i).ne.0) then
      c = c + 1
      candidates(c) = i
    end if
  end do

  ! select a random layer among the occupied ones
  call random_number(die)
  die = 1 + nint(real(c - 1, 8) * die)
  occ_layer = candidates(c)

  return
end function


subroutine random_neighbor(range, x1, x2, L1, L2, n_x1, n_x2)
  implicit none

  integer, intent(in) :: range, x1, x2, L1, L2
  integer, intent(out) :: n_x1, n_x2

  integer :: rmod
  real(kind=8) :: die1, die2

  call random_number(die1)
  call random_number(die2)
  n_x1 = rmod(x1 + (2 * range * nint(die1) - 1), L1)
  n_x2 = rmod(x2 + (2 * range * nint(die2) - 1), L2)
end subroutine


subroutine init_random_seed(m)
  implicit none
  integer, intent(in) :: m
  integer :: i, n
  integer, dimension(:), allocatable :: seed

  call random_seed(size = n)
  allocate(seed(n))

  seed = m + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(PUT = seed)

  deallocate(seed)
end subroutine
