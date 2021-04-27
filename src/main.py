#================================================================================#
#                    Generation of atomistic pseudopotentials                    #
#--------------------------------------------------------------------------------#
#                      Writen from "insert paper reference"                      #
#                               Copyright (c) 2021                               #
#--------------------------------------------------------------------------------#
#                  Author: Ava Dean, University of Oxford, 2021                  #
#================================================================================#

# otfg modules.
import cell
import input_output as io
import generate
import otfg
import parameters
import pspot


currentParams = parameters.get_arguments()

currentCell = cell.UnitCell(currentParams.element)
model = otfg.Model()

soc = False # No spin orbit coupling for now.

# Work out what pseudopotential we are trying to create.
generate.inquire(currentCell, model, soc)

# Update cell ionic charge.
currentCell.ionicCharge = model.ionic_charge

# Initialise pseudopotential.
pseudopotential = pspot.Pseudopotential()

# Allocate variables needed for pseudopotential.
pseudopotential.allocate(currentParams, model)

# Generate the pseudopotential.
generate.generate_psp(currentParams, currentCell, model, pseudopotential)

io.run_time_out(currentParams)
