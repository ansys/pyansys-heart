"""Defines LS-DYNA decks for heart modeling."""

from ansys.dyna.keywords import Deck


class BaseDecks:
    """Class where each attribute corresponds to its respective deck.

    Note
    ----
    Used to distinguish between each of the decks.
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
    """Useful decks for a mechanics simulation."""

    def __init__(self) -> None:
        super().__init__()
        self.cap_elements = Deck()
        self.control_volume = Deck()
        self.pericardium = Deck()


class FiberGenerationDecks(BaseDecks):
    """Useful decks for fiber generation."""

    def __init__(self) -> None:
        super().__init__()
        self.ep_settings = Deck()
        self.create_fiber = Deck()


class PurkinjeGenerationDecks(BaseDecks):
    """Useful decks for Purkinje generation."""

    def __init__(self) -> None:
        super().__init__()
        self.main = Deck()
        self.ep_settings = Deck()


class ElectrophysiologyDecks(BaseDecks):
    """Useful decks for Electrophysiology simulations."""

    def __init__(self) -> None:
        super().__init__()
        self.cell_models = Deck()
        self.ep_settings = Deck()
        self.beam_networks = Deck()


class ElectroMechanicsDecks(ElectrophysiologyDecks, MechanicsDecks):
    """Useful decks for a electromechanics simulation."""

    def __init__(self) -> None:
        super().__init__()
        self.duplicate_nodes = Deck()
