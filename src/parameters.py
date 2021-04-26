import argparse

import data
import input_output as io


def get_arguments():
    parser = argparse.ArgumentParser(prog='otfg',
                                     description='On the fly generation of pseudopotentials')

    parser.add_argument('-e', '--element'   , action='store')
    parser.add_argument('-f', '--functional', action='store')
    parser.add_argument('-s', '--solver'    , action='store')

    args = parser.parse_args()

    return Params(args)


class Params:
    element = None
    energy_unit = 'eV'
    xc_definition = 'PBE'
    relativistic_treatment = 'sh'
    SOC = False

    def __init__(self, args):
        if args.element:
            if args.element.lower() in [e.lower() for e in data.elements]:
                self.element = args.element
            else:
                io.abort('Element {} not recognised.'.format(args.element))
        else:
            io.abort('Please enter an element.')

        if args.functional:
            if args.functional.upper() == 'lda':
                self.xc_definition = 'LDA'
            elif args.functional.upper() == 'PBE':
                self.xc_definition = 'PBE'
            elif args.functional.upper() == 'BLYP':
                self.xc_definition = 'BLYP'
            else:
                io.abort('XC functional {} not known/implemented.'.format(args.functional))
        else:
            self.xc_definition = 'PBE'

        if args.solver:
            if args.solver in ['sh', 'schroedinger']:
                self.relativistic_treatment = 'sh'
            else:
                io.abort('Relativistic treatment {} not known/implemented.'.format(args.solver))
        else:
            self.relativistic_treatment = 'sh'

        self.SOC = True if self.relativistic_treatment == 'dr' else self.SOC
