
class UnitCell:
    def __init__(self, symbol):
        self.speciesSymbol = symbol[0].upper() + symbol[1:].lower()
        self.speciesPot = None
        self.ionicCharge = 0
