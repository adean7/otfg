

class Element:
    def __init__(self, symbol, Z):
        self.symbol = symbol[0].upper() + symbol[1:].lower()
        self.Z = Z

    def __eq__(self, other):
        return True if self.symbol == other.symbol else False


elements = ['H', 'He', 'Li']
periodicTable = { symbol : Element(symbol, Z) for Z, symbol in enumerate(elements, 1)}


# Definitions of pseudopotentials.
ps_def = {"H1"   : "1|0.6|13|15|17|10(qc=8)",       # Hydrogen with pseudopotential 1 (1 corresponds to C19)
          "He1"  : "1|1.0|12|14|16|10(qc=7)",       # Helium with pseudopotential 1 (1 corresponds to C19)
          "Li1"  : "1|1.0|14|16|18|10U:20(qc=7)"    # Lithium with pseudopotential 1 (1 corresponds to C19)
          }

# Groundstate occupations.
gs_occ = { 0 : [10,20,21,30,31,32,40,41,42,43,50,51,52,53,54,60,61,62,70,71], # nl quantum numbers
           1 : [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], # Hydrogen occupancies.
           2 : [ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], # Helium occupancies.
           3 : [ 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Lithium occupancies.
          }

gs_config = { 1 : '1s1' # Hydrogen config.
             }