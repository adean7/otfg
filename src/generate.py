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
    ps_set = '1'    # Corresponds to C19 pseudopotentials for now.

    definition = data.ps_def.get(currentCell.speciesSymbol + ps_set, None)
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
        if indx == -1:
            ctemp = temp_def[0:]
        else:
            ctemp = temp_def[0:indx]

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

        for i in range(len(ctemp)):
            if ctemp[i] == 'N': nb += 2 * l + 1 # Count the norm conserving projectors.
            if ctemp[i] == 'U': nb += 2 * l + 1 # Count the ultrasoft projectors.
            if ctemp[i] == 'H': ng += 2 * l + 1 # Count the norm conserving gamma projectors.
            if ctemp[i] == 'G': ng += 2 * l + 1 # Count the ultrasoft gamma projectors.

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
        if indx == -1:
            ctemp = temp_def[0:]
        else:
            ctemp = temp_def[0:indx]

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
    #define(ab, aeat, pseudopotential)

    # Construct the pseudopotential.
    #construct(ab, aeat, pseudopotential)

    # Descreen the pseudopotential.
    #descreen(ab, aeat, pseudopotential)

    # Write out generated potential.
    #output(ab, aeat, pseudopotential)




def define(ab: basis.AllElectronBasis, aeat: scf.AllElectronAtom, pseudopotential: pspot.Pseudopotential):
    pass

def construct(ab: basis.AllElectronBasis, aeat: scf.AllElectronAtom, pseudopotential: pspot.Pseudopotential):
    pass

def descreen(ab: basis.AllElectronBasis, aeat: scf.AllElectronAtom, pseudopotential: pspot.Pseudopotential):
    pass

def output(ab: basis.AllElectronBasis, aeat: scf.AllElectronAtom, pseudopotential: pspot.Pseudopotential):
    pass

