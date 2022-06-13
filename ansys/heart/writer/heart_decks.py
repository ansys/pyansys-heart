"""Defines LS-DYNA decks for heart modeling."""

from ansys.dyna.keywords import Deck

class BaseDecks:
    """Class where each attribute corresponds to its respective deck. Used to the distinguish between each of the decks.
    This base class defines some commonly used decks.
    """

    def __init__(self) -> None:
        self.main = Deck()
        self.parts = Deck()
        self.nodes = Deck()
        self.solid_elements = Deck()
        self.material = Deck()
        self.segment_sets = Deck()
        self.node_sets = Deck()
        self.boundary_conditions = Deck()        

        return

class MechanicsDecks(BaseDecks):
    """This class inherits from the BaseDecks class and defines additional useful decks"""

    def __init__(self) -> None:
        super().__init__()
        self.cap_elements = Deck()
        self.control_volume = Deck()
        self.pericardium = Deck()

class FiberGenerationDecks(BaseDecks):
    """This class inherits from the BaseDecks class and defines additional useful decks for fiber generation"""
    def __init__(self) -> None:
        super().__init__()
        self.ep_settings = Deck()
        self.create_fiber = Deck()

class ElectrophysiologyDecks(BaseDecks):
    """Adds decks specificly for Electrophysiology simulations"""

    def __init__(self) -> None:
        super().__init__()
