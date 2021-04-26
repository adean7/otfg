#================================================================================#
#                    Generation of atomistic pseudopotentials                    #
#--------------------------------------------------------------------------------#
#                      Writen from "insert paper reference"                      #
#                               Copyright (c) 2021                               #
#--------------------------------------------------------------------------------#
#                  Author: Ava Dean, University of Oxford, 2021                  #
#================================================================================#

# Python modules.
import sys

# otfg modules.
import cell
import generate
import otfg
import pspot


try:
    element = sys.argv[1]
except IndexError:
    print('Please enter a system prefix.')
    sys.exit(1)
else:
    element = element.strip().upper()

currentCell = cell.UnitCell(element)
model = otfg.Model()

soc = False # No spin orbit coupling for now.

# Work out what pseudopotential we are trying to create.
generate.inquire(currentCell, model, soc)

# Update cell ionic charge.
currentCell.ionicCharge = model.ionic_charge

# Initialise pseudopotential.
pseudopotential = pspot.Pseudopotential()

# Allocate variables needed for pseudopotential.
pspot.allocate(pseudopotential, model)

# Generate the pseudopotential.
generate.generate_psp(currentCell, model, pseudopotential)
