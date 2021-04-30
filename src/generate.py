import numpy as np

import basis
import cell
import data
import input_output as io
import otfg
import parameters
import pspot
import scf

import re

def inquire(currentParams, currentCell, model):

    element = data.periodicTable.get(currentCell.speciesSymbol, None)
    if element is None:
        io.abort('Input species {} not known'.format(currentCell.speciesSymbol))

    Z = element.Z

    # Initialise the atomic basis.
    ab = basis.AllElectronBasis()
    ab.init_ae_basis(Z)

    # Get the definition of the pseudopotential.
    ps_set = currentParams.ps_set

    definition = data.ps_def.get(ps_set).get(currentCell.speciesSymbol)
    if definition is None:
        io.abort('Error finding pseudopotential definition for {}, set {}.'.format(currentCell.speciesSymbol, ps_set))

    if definition.count('|') < 5:
        io.abort('otfg definition or library not understood {}.'.format(definition.strip()))

    # Get rid of any spaces.
    definition = definition.replace(' ', '')

    currentCell.speciesPot = definition

    # Set the configuration to the groundstate default.
    cfg_occ = [float(num) for num in data.gs_occ.get(Z, None)]

    # Check if there are any pseudisation flags set.
    if '(' in definition:
        if ')' not in definition:
            io.abort('Syntax error for pseudisation flags.')

        definition = re.sub("[(].*[)]", "", definition)

    # Now proceed with the analysis of the string.
    indx = definition.rfind('|')
    definition = definition[indx+1:]
    if len(definition) < 1:
        io.abort('No reference states specified.')

    nref = 1 + definition.count(':')

    # Find out how many beta functions there are.
    temp_def              = definition
    tot_num_projectors    = 0
    tot_num_bl_projectors = 0
    occ                   = 0

    for n in range(nref):
        indx = temp_def.find(':')
        ctemp = temp_def[0:] if indx == -1 else temp_def[0:indx]

        # Find out the angular momentum of this reference state
        l = int(ctemp[1])

        # Find out the reference state label
        nl = int(ctemp[0:2])

        # Count the number of occupancies.
        for i in range(len(cfg_occ)):
            if nl == data.gs_occ.get(0)[i]:
                occ += cfg_occ[i]

        temp_def = temp_def[indx+1:]
        nb = 0
        ng = 0

        for char in ctemp.strip():
            if char == 'N': nb += 2 * l + 1 # Count the norm conserving projectors.
            if char == 'U': nb += 2 * l + 1 # Count the ultrasoft projectors.
            if char == 'H': ng += 2 * l + 1 # Count the norm conserving gamma projectors.
            if char == 'G': ng += 2 * l + 1 # Count the ultrasoft gamma projectors.

        if nb == 0 and ctemp.find('P') == -1 and ctemp.find('L') == -1:
            nb = 2 * (2 * l + 1)

        if currentParams.SOC and l > 0: # Double for j
            tot_num_projectors    += 2 * nb
            tot_num_bl_projectors += 2 * ng
        else:
            tot_num_projectors    += nb
            tot_num_bl_projectors += ng

    tot_num_bl_projectors += tot_num_projectors

    # For ELNES calculation.
    num_core_projectors = 1  # We assume there is at least one core projector. Even for H.
    temp_def = definition
    nnmax = -1
    for n in range(nref):
        indx = temp_def.find(':')
        ctemp = temp_def[0:] if indx == -1 else temp_def[0:indx]

        # Find out the principal quantum number label
        nn = int(ctemp[0])

        if n == 1:
            nnmax = nn
        else:
            nnmax = max(nn, nnmax)

        temp_def = temp_def[indx+1:]

    # We consider configuration of core orbitals is a noble config. .i.e. it is 1s, or 1s,2s,2p,
    # or 1s,2s,2p,3s,3p, or etc. ---ELNES

    nnmax = min(6, nnmax)

    if nnmax in [1, 2]:
        num_core_projectors = 1  # 1s
    elif nnmax == 3:
        num_core_projectors = 5  # 1s, 2s, 2p
    elif nnmax == 4:
        num_core_projectors = 9  # 1s, 2s, 2p, 3s, 3p
    elif nnmax == 5:
        num_core_projectors = 18 # 1s, 2s, 2p, 3s, 3p, 3d, 4s, 4p
    elif nnmax == 6:
        num_core_projectors = 27 # 1s, 2s, 2p, 3s, 3p, 3d, 4s, 4p, 4d, 5s, 5p
    elif nnmax == 7:
        num_core_projectors = 43 # 1s, 2s, 2p, 3s, 3p, 3d, 4s, 4p, 4d, 5s, 5p, 4f, 5d, 6s, 6p
    else:
        io.abort('Error in setting core_projectors: there are no projectors?')

    # Set the other parameters for this species.
    tot_num_points = 2001
    ionic_charge = float(Z) - sum(cfg_occ) + occ
    num_mesh = ab.npts
    num_kkbeta = ab.npts
    num_nqf = 10
    bl_points = tot_num_points # Same as beta projectors

    model.update_info(tot_num_projectors,
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
                      cfg_occ)

    # Update cell ionic charge.
    currentCell.ionicCharge = model.ionic_charge



def generate_psp(currentParams: parameters.Params, currentCell    :  cell.UnitCell,
                 model        :       otfg.Model , pseudopotential: pspot.Pseudopotential):
    Z = model.Z
    definition = currentCell.speciesPot

    # Set the configuration to the groundstate default.
    cfg_occ = model.cfg_occ

    # Store the generation configuration.
    pseudopotential.gen_cfg = cfg_occ
    cfg_occ_save = cfg_occ   # ELNES

    # Get the atomic basis.
    ab = model.ab

    # Initialise the all electron atom.
    aeat = scf.AllElectronAtom()
    aeat.atomic_solver = pseudopotential.atomic_solver
    aeat.SOC = pseudopotential.SOC
    aeat.theory = pseudopotential.theory
    aeat.init_atom(Z, ab, cfg_occ)

    # Set the all electron XC used.
    aeat.theory = pseudopotential.theory

    # Set the finite nucleus.
    aeat.fin_nuc = True if aeat.atomic_solver in ['kh', 'sh'] else False

    # Do the all electron calculation for the reference atom.
    scf.cycle(currentParams, currentCell, ab, aeat)

    # Make the ZORA K function.
    # Don't need this just yet as only have Schroedinger solver for now.

    # Define the pseudopotential to be generated.
    #define(currentParams, currentCell, ab, aeat, pseudopotential)

    # Construct the pseudopotential.
    #construct(ab, aeat, pseudopotential)

    # Descreen the pseudopotential.
    #descreen(ab, aeat, pseudopotential)

    # Write out generated potential.
    #output(ab, aeat, pseudopotential)




def define(currentParams: parameters.Params,
           currentCell: cell.UnitCell,
           ab: basis.AllElectronBasis,
           aeat: scf.AllElectronAtom,
           ps: pspot.Pseudopotential):

    # Define comment characters.
    comment = ['#', ';', '!', '^M']

    # Set the working pseudopotential definition.
    definition = currentCell.speciesPot

    # We have the pseudopotential set in pseudopotential.set_

    # Output testing files? - not doing this anyway.
    testing = False if definition.find('[') == -1 else True

    # Check if there are any pseudisation flags set.
    if '(' in definition:
        if ')' not in definition:
            io.abort('Syntax error for pseudisation flags.')

        ps.flags = definition[definition.find('(') + 1:definition.find(')')].strip()
        definition = re.sub("[(].*[)]", "", definition)

    if '{' in definition:
        if '}' not in definition:
            io.abort('Syntax error for config adjust.')
        definition = re.sub("[{].*[}]", "", definition)

    if '[' in definition:
        if ']' not in definition:
            io.abort('Syntax error for config adjust.')
        definition = re.sub("[\[].*[]]", "", definition)

    # Set the local channel.
    ps.pv = False      # Don't do the polynomial fit by default.

    indx = definition.find('|')
    temp = definition[:indx]

    if temp.find('=') != -1:
        io.abort('Overide mode not implemented.')
    else:
        overide = False
        ps.local_e = 0.0   # The energy.

        try:
            ps.local_l = abs(int(temp))
        except Exception:
            io.abort('Problem parsing definition 01.')

        if temp.find('-') != -1:
            ps.pv = True

    ps.local_norm = 0      # Do not conserve the norm.

    # Find out about the core radii - if only one is specified, use defaults
    # - note, they may be specified on a projector by projector basis in addition.

    # Count the number of "|" -- if 5, then use defaults.
    n = definition.count('|')

    # Read the core radii in from the definition string.
    jndx = definition.find('|') + 1

    if n == 5:      # Only one rc specified - use defaults.
        indx = definition[jndx:].find('|')
        try:
            ps.local_rc = float(definition[jndx:jndx + indx])
        except Exception:
            io.abort('Problem parsing definition 02.')

        temp_rc = ps.local_rc
        ps.rinner = temp_rc * 0.7      # Default setting.
        definition = definition[jndx+1+indx:]       # +1 to get rid of next |

    elif n == 7:    # local_rc, beta_rc
        indx = definition[jndx:].find('|')
        try:
            ps.local_rc = float(definition[jndx:jndx + indx])
        except Exception:
            io.abort('Problem parsing definition 03.')

        jndx += 1 + indx

        indx = definition[jndx:].find('|')
        try:
            temp_rc = float(definition[jndx:jndx+indx])
        except Exception:
            io.abort('Problem parsing definition 04.')

        jndx += 1 + indx

        indx = definition[jndx:].find('|')
        try:
            ps.rinner[0] = float(definition[jndx:jndx + indx])
        except Exception:
            io.abort('Problem parsing definition 05.')
        for i in range(len(ps.rinner)):
            ps.rinner[i] = ps.rinner[0]

        definition = definition[jndx+1+indx:]  # +1 to get rid of next |

    else:
        io.abort('Unable to parse pseudopotential - incorrect number of "|".')

    temp_rc2 = temp_rc      # Keep a copy of default beta_rc to use for the gamma projectors.

    # Read in the definitions of COARSE, MEDIUM and FINE
    # - note, the should be in atomic units (Ha).
    indx = definition.find('|')

    try:
        ps.coarse = float(definition[:indx])
    except:
        io.abort('Problem parsing definition 06.')

    jndx = 1 + indx
    indx = definition[jndx:].find('|')

    try:
        ps.medium = float(definition[jndx:jndx + indx])
    except:
        io.abort('Problem parsing definition 07.')

    jndx += 1 + indx
    indx = definition[jndx:].find('|')

    try:
        ps.fine = float(definition[jndx:jndx + indx])
    except:
        io.abort('Problem parsing definition 08.')

    definition = definition[jndx+1+indx:]  # +1 to get rid of next |

    # Set PRECISE and EXTREME (we don't really know what EXTREME is).
    ps.precise = max(1.03 * ps.fine + 1.57, 1.2 * ps.fine)
    ps.extreme = 1.6 * ps.fine

    # Set the cutoff energy.
    if currentParams.basis_precision == 'COARSE':
        currentParams.cut_off_energy = max(ps.coarse, currentParams.cut_off_energy)
    elif currentParams.basis_precision == 'MEDIUM':
        currentParams.cut_off_energy = max(ps.medium, currentParams.cut_off_energy)
    elif currentParams.basis_precision == 'FINE':
        currentParams.cut_off_energy = max(ps.fine, currentParams.cut_off_energy)
    elif currentParams.basis_precision == 'PRECISE':
        currentParams.cut_off_energy = max(ps.precise, currentParams.cut_off_energy)
    elif currentParams.basis_precision == 'EXTREME':
        currentParams.cut_off_energy = max(ps.extreme, currentParams.cut_off_energy)
    else:
        io.abort('Don\'t recognise basis precision.')

    # Now read in the projector definitions.

    # Count the number of reference states.
    nref = definition.count(':') + 1

    # Find out how many beta functions there are so the arrays can be allocated.
    temp_def = definition
    ps.num_beta = 0
    for n in range(nref):
        indx = temp_def.find(':')
        ctemp = temp_def[:indx] if indx > -1 else temp_def
        ctemp = ctemp.strip()
        temp_def = temp_def[indx+1:]

        nb = 0
        for char in ctemp:
            if char == 'N': nb += 1     # Count the norm conserving projectors.
            if char == 'U': nb += 1     # Count the ultrasoft projectors.

        nb = 2 if nb == 0 else nb

        for char in ctemp:
            if char == 'L': nb = 0      # Don't count the local channel.
            if char == 'P': nb = 0      # Don't count the dummy channel.

        ps.num_beta += nb

    # Allocate the beta arrays.
    ps.beta_l = np.zeros(ps.num_beta, dtype=int)
    ps.beta_e = np.zeros(ps.num_beta)
    ps.beta_rc = np.zeros(ps.num_beta)
    ps.beta_delta = np.zeros(ps.num_beta)
    ps.beta_nrc = np.zeros(ps.num_beta, dtype=int)
    ps.beta_norm = np.zeros(ps.num_beta, dtype=int)
    ps.beta_scheme = np.zeros(ps.num_beta, dtype=str)

    # Set the beta arrays.
    temp_def = definition
    nbb = 0

    for n in range(nref):
        indx = temp_def.find(':')
        ctemp = temp_def[:indx] if indx > -1 else temp_def
        ctemp = ctemp.strip()
        temp_def = temp_def[indx+1:]

        try:
            no, l = int(ctemp[0]), int(ctemp[1])
        except Exception:
            print('Problem parsing definition 09.')

        ps.pseudise[no-1-l,l] = True

        try:
            if ctemp[2] == 'P':
                continue        # Just pseudise - no projectors etc.
        except IndexError:
            pass

        if len(ctemp) > 2 and ctemp[2] == 'L':
            if l != ps.local_l:
                io.abort('Local channel conflict.')

            if not overide:
                ps.local_e = aeat.eval[no-1-l,l]

            ps.local_norm = 0 if ctemp[3] == 'G' else 1

        else:
            ps.beta_e[nbb+1] = aeat.eval[no-1-l,l]

            nb = 0
            for i, char_i in enumerate(ctemp):
                if char_i == 'N' or char_i == 'U':
                    nb += 1     # Count the projectors.

                    for j, char_j in enumerate(ctemp[i:]):
                        if char_j in ['N', 'U', 'H', 'G']:
                            break

                    print(i, j, ctemp)
                    ctemp2 = ctemp[i:j]
                    print(ctemp2)
                    exit()

























































def construct(ab: basis.AllElectronBasis, aeat: scf.AllElectronAtom, pseudopotential: pspot.Pseudopotential):
    pass

def descreen(ab: basis.AllElectronBasis, aeat: scf.AllElectronAtom, pseudopotential: pspot.Pseudopotential):
    pass

def output(ab: basis.AllElectronBasis, aeat: scf.AllElectronAtom, pseudopotential: pspot.Pseudopotential):
    pass

