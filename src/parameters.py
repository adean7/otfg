import argparse
import time

import data
import input_output as io


def get_arguments():
    parser = argparse.ArgumentParser(prog='otfg',
                                     description='On the fly generation of pseudopotentials')

    parser.add_argument('-e', '--element'   , action='store')
    parser.add_argument('-f', '--functional', action='store')
    parser.add_argument('-s', '--string',     action='store')
    parser.add_argument('-S', '--solver'    , action='store')

    args = parser.parse_args()

    return Params(args)


class Params:
    # Set defaults here.
    element                = None
    energy_unit            = 'eV'
    xc_definition          = 'PBE'
    ps_set                 = 'C19'
    basis_precision        = 'FINE'
    cut_off_energy         = 50.0
    relativistic_treatment = 'sh'
    SOC                    = False

    start_time             = time.time()

    def __init__(self, args):
        self.set_element(args)
        self.set_functional(args)
        self.set_string(args)
        self.set_solver(args)
        self.set_SOC()

    def set_element(self, args):
        if args.element:
            if args.element.upper() in [e.upper() for e in data.elements]:
                self.element = args.element[0].upper() + args.element[1:].lower()
            elif args.element in [str(num) for num in range(1, len(data.elements)+1)]:
                el_str = data.elements[int(args.element)-1]
                self.element = el_str[0].upper() + el_str[1:].lower()
            else:
                io.abort('Element {} not recognised.'.format(args.element))
        else:
            io.abort('Please enter an element.')

    def set_functional(self, args):
        if args.functional:
            if args.functional.upper() == 'LDA':
                self.xc_definition = 'LDA'
            elif args.functional.upper() == 'PBE':
                self.xc_definition = 'PBE'
            elif args.functional.upper() == 'BLYP':
                self.xc_definition = 'BLYP'
            else:
                io.abort('XC functional {} not known/implemented.'.format(args.functional))

    def set_string(self, args):
        if args.string:
            if args.string.upper() == 'C19':
                self.ps_set = 'C19'
            else:
                io.abort('Pseudopotential set {} not known/implemented.'.format(args.string))

    def set_solver(self, args):
        if args.solver:
            if args.solver.upper() in ['SH', 'SCHROEDINGER']:
                self.relativistic_treatment = 'sh'
            else:
                io.abort('Relativistic treatment {} not known/implemented.'.format(args.solver))

    def set_SOC(self):
        self.SOC = True if self.relativistic_treatment == 'dr' else self.SOC

        # No SOC for now.
        if self.SOC:
            print('Warning: setting SOC to False as not implemented yet.')
            self.SOC = False

    def get_run_time(self):
        return time.time() - self.start_time