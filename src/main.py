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


# Get the parameters for the calculation.
currentParams = parameters.get_arguments()

# Get the current unit cell.
currentCell = cell.UnitCell(currentParams.element)

# Create a class to hold all the information
model = otfg.Model()

# Work out what pseudopotential we are trying to create.
generate.inquire(currentParams, currentCell, model)

# Initialise pseudopotential.
pseudopotential = pspot.Pseudopotential()

# Allocate variables needed for pseudopotential.
pseudopotential.allocate(currentParams, model)

# Generate the pseudopotential.
generate.generate_psp(currentParams, currentCell, model, pseudopotential)

io.run_time(currentParams)
