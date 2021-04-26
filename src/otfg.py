

class Model:
    num_projectors     : int   = None
    num_core_projectors: int   = None
    num_bl_projectors  : int   = None
    num_points         : int   = None
    num_bl_points      : int   = None
    num_mesh           : int   = None
    num_kkbeta         : int   = None
    num_nqf            : int   = None
    ionic_charge       : float = None
    Z                  : int   = None
    ab                         = None
    cfg_occ                    = None

    def update_info(self,
                    tot_num_projectors,
                    num_core_projectors,
                    tot_num_bl_projectors,
                    tot_num_points,
                    bl_points,
                    num_mesh,
                    num_kkbeta,
                    num_nqf,
                    ionic_charge,
                    Z,
                    ab,
                    cfg_occ):

        self.num_projectors      = tot_num_projectors
        self.num_core_projectors = num_core_projectors
        self.num_bl_projectors   = tot_num_bl_projectors
        self.num_points          = tot_num_points
        self.num_bl_points       = bl_points
        self.num_mesh            = num_mesh
        self.num_kkbeta          = num_kkbeta
        self.num_nqf             = num_nqf
        self.ionic_charge        = ionic_charge
        self.Z                   = Z
        self.ab                  = ab
        self.cfg_occ             = cfg_occ