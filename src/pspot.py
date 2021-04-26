import numpy as np

def allocate(currentParams, pseudopotential, model):
    lmax = 4

    pseudopotential.ID = 0
    pseudopotential.pseudise = np.zeros((10, 5))

    # The ionic charge of the pseudopotential.
    pseudopotential.ionic_charge = 0.0

    # The generation configuration.
    pseudopotential.gen_cfg = np.zeros(20)

    # Number of points in reciprocal radial grid.
    pseudopotential.num_points = 0

    # Maximum g-vector in reciprocal radial grid.
    pseudopotential.gmax = 0.0

    # Local part of ionic pseudopotential.
    pseudopotential.V_loc = np.zeros(model.num_points)

    # Local part of ionic psp, in real space x r.
    pseudopotential.V_locr = np.zeros(model.num_mesh)

    # Derivative of local part of ionic psp wrt g.
    pseudopotential.dV_loc_dg = np.zeros(model.num_points)

    pseudopotential.has_V_tau = False

    # Local part of ionic MGGA in reciprocal space.
    pseudopotential.V_tau = np.zeros(model.num_points)

    # Local part of ionic MGGA in real space.
    pseudopotential.V_taur = np.zeros(model.num_mesh)

    # Derivative of local part of ionic MGGA wrt g.
    pseudopotential.dV_tau_dg = np.zeros(model.num_points)

    # Number of projectors.
    pseudopotential.num_ps_projectors = 0

    # q_nm for each pair of projectors.
    pseudopotential.q = np.zeros((model.num_projectors, model.num_projectors))

    # D^0_nm for each pair of projectors.
    pseudopotential.D0 = np.zeros((model.num_projectors, model.num_projectors))

    # Assume the pseudopotential does not include spin orbit coupling unless otherwise specified.
    pseudopotential.SOC = currentParams.SOC

    # Angular momentum, l, for each projector.
    pseudopotential.projector_l = np.zeros(model.num_projectors, dtype=int)

    # Angular momentum, m, for each projector.
    pseudopotential.projector_m = np.zeros(model.num_projectors, dtype=int)

    # Twice the total angular momentum, j, for each projector.
    pseudopotential.projector_2j = np.zeros(model.num_projectors, dtype=int)

    # Is there a core charge for calculating NLCC?.
    pseudopotential.core_charge = True

    # The radial core charge density.
    pseudopotential.rs_core_charge = np.zeros(model.num_mesh)
    pseudopotential.core_radial_charge = np.zeros(model.num_points)

    # Derivative of the radial core charge density.
    pseudopotential.dcore_radial_charge_dg = np.zeros(model.num_points)

    # The radial core kinetic_energy_density.
    pseudopotential.rs_core_kinetic_energy_density = np.zeros(model.num_mesh)
    pseudopotential.core_radial_kinetic_energy_density = np.zeros(model.num_points)

    # Derivative of the radial core kinetic energy density.
    pseudopotential.dcore_radial_kinetic_energy_density_dg = np.zeros(model.num_points)

    # Non-local part of ionic psp in recip space.
    pseudopotential.beta = np.zeros((model.num_points, model.num_projectors))

    # Non-local part of ionic psp on a real log grid.
    pseudopotential.betar = np.zeros((model.num_mesh, model.num_projectors))

    # Non-local part of ionic psp in real space.
    pseudopotential.beta_real = np.zeros((model.num_points, model.num_projectors))

    # Magnitude of maximum radius of radial grid.
    pseudopotential.rmax = 0.0

    # The real space radii.
    pseudopotential.real_radii = 0.0

    # The size of the radial grid for the aug fns.
    pseudopotential.mesh = 0

    # The extent of the aug fns.
    pseudopotential.kkbeta = 0

    # The number of optimisation functions.
    pseudopotential.nqf = 0

    # The radial logarithmic grid.
    pseudopotential.rlog = np.zeros(model.num_mesh)
    pseudopotential.rab = np.zeros(model.num_mesh)

    # The inner core radii for the aug fns.
    pseudopotential.rinner = np.zeros(2 * lmax + 1)
    pseudopotential.ninner = np.zeros(2 * lmax + 1)

    # The aug fns on a radial real space grid.
    pseudopotential.qfunc = np.zeros((model.num_kkbeta, model.num_projectors, model.num_projectors))

    # The kinetic energy density aug fns on a radial real space grid.
    pseudopotential.tfunc = np.zeros((model.num_kkbeta, 2 * lmax + 1, model.num_projectors, model.num_projectors))

    # The polynomial coefficients of the aug fns.
    pseudopotential.qfcoef = np.zeros((model.num_nqf, 2 * lmax + 1, model.num_projectors, model.num_projectors))

    # The polynomial coefficients of the aug fns of the kinetic energy density.
    pseudopotential.tfcoef = np.zeros((model.num_nqf, 2 * lmax + 1, model.num_projectors, model.num_projectors))

    # The l-dependent augmentation functions on a radial grid.
    pseudopotential.qfunc_with_l = np.zeros((model.num_kkbeta, 2 * lmax + 1, model.num_projectors, model.num_projectors))

    pseudopotential.qfunc_is_polynomial = True # pspot defaults to augmentation functions being stored as polynomial expansions

    # The extra information stored as flags.
    pseudopotential.flags = ''

    #--------------------------------------------------
    # Data associated with the Blochl gamma projectors.
    #--------------------------------------------------

    # Non-local part of ionic psp in recip space
    pseudopotential.gamma = np.zeros((model.num_points, model.num_bl_projectors))

    # Non-local part of ionic psp in real space.
    # Angular momentum, l, for each blochl projector.
    pseudopotential.bl_projector_l = np.zeros(model.num_bl_projectors)

    # Angular momentum, m, for each blochl projector.
    pseudopotential.bl_projector_m = np.zeros(model.num_bl_projectors)

    # Number of gamma projectors.
    pseudopotential.num_bl_projectors = 0

    # Number of points in reciprocal radial grid.
    pseudopotential.bl_num_points = 0

    # Maximum g-vector in reciprocal radial grid.
    pseudopotential.bl_gmax = 0.0

    #--------------------------------------------------
    # Data associated with the NMR Blochl Augmentation.
    #--------------------------------------------------

    # Initialise the core shift and susceptibility.
    pseudopotential.core_nmr_shift = 0.0
    pseudopotential.core_mag_sus = 0.0

    # The augmentation factors for each projector pair.
    pseudopotential.bl_dia_nmrshift = np.zeros((model.num_bl_projectors, model.num_bl_projectors))
    pseudopotential.bl_para_nmrshift = np.zeros((model.num_bl_projectors, model.num_bl_projectors))
    pseudopotential.bl_para_damped = np.zeros((model.num_bl_projectors, model.num_bl_projectors))
    pseudopotential.bl_efg_aug = np.zeros((model.num_bl_projectors, model.num_bl_projectors))
    pseudopotential.bl_hf_fermi = np.zeros((model.num_bl_projectors, model.num_bl_projectors))
    pseudopotential.bl_hf_dipolar = np.zeros((model.num_bl_projectors, model.num_bl_projectors))
    pseudopotential.bl_hf_dipolar_cutoff = np.zeros((model.num_bl_projectors, model.num_bl_projectors))
    pseudopotential.bl_hf_zora = np.zeros((model.num_bl_projectors, model.num_bl_projectors))

    # Data associated with EPR g-tensor.
    pseudopotential.bl_zke_aug_gtensor = np.zeros((model.num_bl_projectors, model.num_bl_projectors))
    pseudopotential.bl_sodia_g_tensor = np.zeros((model.num_bl_projectors, model.num_bl_projectors))
    pseudopotential.bl_sopara_g_tensor = np.zeros((model.num_bl_projectors, model.num_bl_projectors))

    #--------------------------------------------------------------------------
    # Data associated with augmentation terms for core excited matrix elements.
    #--------------------------------------------------------------------------

    # Number of core projectors.
    pseudopotential.num_core_projectors = 0

    # Number of core radial orbitals.
    pseudopotential.num_core_radials = 0

    # Principal quantum number of core radials.
    pseudopotential.core_radials_n = np.zeros(15, dtype=int)

    # Angular momentum l, for each core radials.
    pseudopotential.core_radials_l = np.zeros(15, dtype=int)

    # The radial core orbitals (reciprocal space).
    pseudopotential.core_radials = np.zeros((model.num_points, 15, 5))

    # The ELNES augmentation term.
    pseudopotential.bl_elnes_aug = np.zeros((model.num_core_projectors, model.num_bl_projectors))

    # Just LDA and PBE for now.
    pseudopotential.theory = currentParams.xc_definition

    # Just Schroedinger for now.
    pseudopotential.atomic_solver = currentParams.relativistic_treatment






class Pseudopotential:
    Z                                      = None #integer                          # The atomic number of the atom
    ionic_charge                           = None #real(kind=dp)                    # The ionic charge of the pseudopotential
    gen_cfg                                = None #real(kind=dp), dimension(20)     # The configuration used in generation
    set_                                   = None #character(len=20)                # The pseudopotential set used
    theory                                 = None #character(len=20)                # The level of theory used
    atomic_solver                          = None #character(len=2)                 # The atomic solver used (Schroedinger , KH etc)
    pseudise                               = None #logical, dimension(:,:)          # Which channels to be pseudised
    num_points                             = None #integer                          # Number of points in reciprocal radial grid
    gmax                                   = None #real(kind=dp)                    # Maximum g-vector in reciprocal radial grid
    pv                                     = None #logical                          # Do the polynomial fit to the potential
    local_l                                = None #integer                          # Angular momentum channel of V_loc
    local_e                                = None #real(kind=dp)                    # The energy for the local partial wave
    local_rc                               = None #real(kind=dp)                    # Cutoff for the local potential
    local_norm                             = None #integer                          # Apply norm contraints?
    local_scheme                           = None #character(len=2)                 # Pseudisation scheme used
    local_nrc                              = None #integer                          # Cutoff point for the local potential
    V_loc                                  = None #real(kind=dp), dimension(:)      # Local part of ionic psp
    V_locr                                 = None #real(kind=dp), dimension(:)      # Local part of ionic psp, in real space x r
    has_V_tau                              = None #logical                          # Flag to determine presence of V_tau
    V_tau                                  = None #real(kind=dp), dimension(:)
    V_taur                                 = None #real(kind=dp), dimension(:)
    dV_loc_dg                              = None #real(kind=dp), dimension(:)      # Derivative of local part of ionic psp wrt g
    dV_tau_dg                              = None #real(kind=dp), dimension(:)      # Derivative of local part of ionic MGGA wrt g
    num_beta                               = None #integer                          # Number of projectors (l only)
    num_ps_projectors                      = None #integer                          # Number of projectors (l and m)
    q                                      = None #real(kind=dp), dimension(:,:)    # q_nm for each pair of projectors
    D0                                     = None #real(kind=dp), dimension(:,:)    # D^0_nm for each pair of projectors
    SOC                                    = None #logical                          # Whether or not this pseudopotential includes spin orbit coupling.
    beta_l                                 = None #integer,       dimension(:)      # Angular momentum for each projector (l only)
    beta_2j                                = None #integer,       dimension(:)      # Total angular momentum for each projector (l only)
    beta_e                                 = None #real(kind=dp), dimension(:)      # The energies for the projector partial waves
    beta_rc                                = None #real(kind=dp), dimension(:)      # Cutoff radii for pseudising the projectors
    beta_delta                             = None #real(kind=dp), dimension(:)      # The eigenvalue shift for this channel
    beta_nrc                               = None #integer,       dimension(:)      # Cutoff points for pseudising the projectors
    beta_norm                              = None #integer,       dimension(:)      # Apply norm contraints to the projectors?
    beta_scheme                            = None #character(len=2),dimension(:)    # Pseudisation scheme used
    projector_l                            = None #integer,       dimension(:)      # Angular momentum, l, for each projector
    projector_m                            = None #integer,       dimension(:)      # Angular momentum, m, for each projector
    projector_2j                           = None #integer,       dimension(:)      # Twice the total angular momentum, j, for each projector
    core_charge                            = None #logical                          # Is there a core charge for calculating NLCC?
    rs_core_charge                         = None #real(kind=dp), dimension(:)      # The radial core charge density in realspace
    core_radial_charge                     = None #real(kind=dp), dimension(:)      # The radial core charge density
    dcore_radial_charge_dg                 = None #real(kind=dp), dimension(:)      # Derivative of the radial core charge density
    rs_core_kinetic_energy_density         = None #real(kind=dp), dimension(:)      # The radial core kinetic energy density in realspace
    core_radial_kinetic_energy_density     = None #real(kind=dp), dimension(:)      # The radial core kinetic energy density
    dcore_radial_kinetic_energy_density_dg = None #real(kind=dp), dimension(:)      # Derivative of the radial core kinetic energy density
    beta                                   = None #real(kind=dp), dimension(:,:)    # Non-local part of ionic psp in recip space
    betar                                  = None #real(kind=dp), dimension(:,:)    # Non-local part of ionic psp on log real grid
    beta_real                              = None #real(kind=dp), dimension(:,:)    # Non-local part of ionic psp in real space
    rmax                                   = None #real(kind=dp)                    # Magnitude of maximum radius of radial grid
    rcmax                                  = None #real(kind=dp)                    # Magnitude of the high core radius
    real_radii                             = None #real(kind=dp)                    # The real space radii
    mesh                                   = None #integer                          # The size of the radial grid for the aug fns
    kkbeta                                 = None #integer                          # The extent of the aug fns
    nqpts                                  = None #integer                          # The number of points for aug fn optimisation
    qc                                     = None #real(kind=dp)                    # The optimisation q cutoff
    nqf                                    = None #integer                          # The number of optimisation functions
    nderiv                                 = None #integer                          # The number of derivatives to match
    nextra                                 = None #integer                          # The number of extra degrees of freedom
    rlog                                   = None #real(kind=dp), dimension(:)      # The radial logarithmic grid
    rab                                    = None #real(kind=dp), dimension(:)      # The radial logarithmic grid
    rinner                                 = None #real(kind=dp), dimension(:)      # The inner core radii for the aug fns
    ninner                                 = None #integer,       dimension(:)      # The inner core cutoff points for the aug fns
    rnlcc                                  = None #real(kind=dp)                    # The core radii for the core charge
    nnlcc                                  = None #integer                          # The core cutoff points for the core charge
    qfunc                                  = None #real(kind=dp), dimension(:,:,:)  # The aug fns on a radial real space grid
    tfunc                                  = None #real(kind=dp), dimension(:,:,:,:)# The kinetic energy density aug fns on a radial real space grid
    qfcoef                                 = None #real(kind=dp), dimension(:,:,:,:)# The polynomial coefficients of the aug fns
    tfcoef                                 = None #real(kind=dp), dimension(:,:,:,:)# The polynomial coefficients of the aug fns
    qfunc_with_l                           = None #real(kind=dp), dimension(:,:,:,:)# The l-dependent augmentation functions on a radial grid
    qfunc_is_polynomial                    = None #logical                          # Are augmentation functions for r<rinner saved as polynomial coeffs
    flags                                  = None #character(len=40)                # Flags for adjusting the pseudisation scheme
    coarse                                 = None #real(kind=dp)                    # The recommended coarse cutoff for this species
    medium                                 = None #real(kind=dp)                    # The recommended medium cutoff for this species
    fine                                   = None #real(kind=dp)                    # The recommended fine cutoff for this species
    precise                                = None #real(kind=dp)                    # The recommended precise cutoff for this species
    extreme                                = None #real(kind=dp)                    # The recommended extreme cutoff for this species
    calc_c                                 = None #real(kind=dp)                    # The derived coarse cutoff for this species
    calc_m                                 = None #real(kind=dp)                    # The derived medium cutoff for this species
    calc_f                                 = None #real(kind=dp)                    # The derived fine cutoff for this species
    calc_e                                 = None #real(kind=dp)                    # The derived fine cutoff for this species
    calc_p                                 = None #real(kind=dp)                    # The derived precise cutoff for this species
    aug                                    = None #logical                          # Do any augmentations?
    scpsp                                  = None #logical                          # Is this an self consistent pseudopotential?
    ID                                     = None #integer                          # An id number for the pseudopotential

    #-------------------------------------------------
    # Data associated with the Blochl gamma projectors
    #-------------------------------------------------

    gamma_l                                = None #integer,       dimension(:)      # Angular momentum for each projector (l only)
    gamma_2j                               = None #integer,       dimension(:)      # Total angular momentum for each projector
    gamma_e                                = None #real(kind=dp), dimension(:)      # The energies for the projector partial waves
    gamma_rc                               = None #real(kind=dp), dimension(:)      # Cutoff radii for pseudising the projectors
    gamma_scheme                           = None #character(len=2), dimension(:)   # Pseudisation scheme used
    gamma_nrc                              = None #integer,       dimension(:)      # Cutoff points for pseudising the projectors
    gamma_norm                             = None #integer,       dimension(:)      # Apply norm contraints to the projectors?
    num_gamma                              = None #integer                          # Number of projectors (l only)
    num_bl_projectors                      = None #integer                          # Number of projectors
    gamma                                  = None #real(kind=dp), dimension(:,:)    # Non-local part of ionic psp in recip space
    gamma_real                             = None #real(kind=dp), dimension(:,:)    # Non-local part of ionic psp in real space
    bl_projector_l                         = None #integer,       dimension(:)      # Angular mom, l, for each projector
    bl_projector_m                         = None #integer,       dimension(:)      # Angular momentum, m, for each projector
    max_real_bl_points                     = None #integer                          # Max number of real space pts for a given gamma
    bl_gmax                                = None #real(kind=dp)                    # Magnitude of maximum g-vector in radial grid
    bl_num_points                          = None #integer                          # Num pts in recip radial grid for each species

    #-------------------------------------------------
    # Data associated with the NMR Blochl Augmentation
    #-------------------------------------------------

    core_nmr_shift                         = None #real(kind=dp)                    # The core chemical shielding for each species
    core_mag_sus                           = None #real(kind=dp)                    # The core magnetic susceptibility for each species
    bl_dia_nmrshift                        = None #real(kind=dp), dimension(:,:)    # The diamagnetic augmentation
    bl_para_nmrshift                       = None #real(kind=dp), dimension(:,:)    # The paramagnetic augmentation
    bl_para_damped                         = None #real(kind=dp), dimension(:,:)    # The paramagnetic augmentation, with damping for JC
    bl_efg_aug                             = None #real(kind=dp), dimension(:,:)    # The electric field gradient augmentation
    bl_hf_fermi                            = None #real(kind=dp), dimension(:,:)    # The hf coupling fermi contact augmentation
    bl_hf_dipolar                          = None #real(kind=dp), dimension(:,:)    # The hf coupling dipolar augmentation
    bl_hf_dipolar_cutoff                   = None #real(kind=dp), dimension(:,:)    # The hf coupling dipolar augmentation when the dipolar term is damped
    bl_hf_zora                             = None #real(kind=dp), dimension(:,:)    # The hf coupling ZORA augmentation (ZSO)

    #----------------------------------
    # Data associated with EPR g-tensor
    #----------------------------------

    ke_aug_gtensor                         = None #real(kind=dp), dimension(:,:)    # The rmc augmentation factor for each projector pair
    dia_so_gtensor                         = None #!real(kind=dp), dimension(:,:)   # The diamagnetic so g_tensor factor
    odia_g_tensor                          = None #real(kind=dp), dimension(:,:)    # The spin orbit diamagnetic augmentation factor for each projector pair
    opara_g_tensor                         = None #real(kind=dp), dimension(:,:)    # The spin orbit paramagnetic augmentation factor for each projector pair

    #--------------------------------------------------
    # Data associated with augmentation terms for ELNES
    #--------------------------------------------------

    num_core_projectors                    = None #integer                          # Number of core projectors
    num_core_radials                       = None #integer                          # Number of core radial orbitals
    core_radials_n                         = None #integer,       dimension(:)      # Principal quantum number of core radials
    core_radials_l                         = None #integer,       dimension(:)      # Angular momentum l, for each core radials
    core_radials                           = None #real(kind=dp), dimension(:,:,:)  # The radial core orbitals (reciprocal space)
    bl_elnes_aug                           = None #real(kind=dp), dimension(:,:)    # The elnes augmentation
