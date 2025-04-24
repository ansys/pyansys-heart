# Copyright (C) 2023 - 2025 ANSYS, Inc. and/or its affiliates.
# SPDX-License-Identifier: MIT
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Defines LS-DYNA decks for heart modeling."""

from ansys.dyna.core import Deck


class BaseDecks:
    """Class where each attribute corresponds to its respective deck.

    Notes
    -----
    This class used to distinguish between each of the decks.
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

    def add_deck(self, deckname: str):
        """Add deck by filename."""
        setattr(self, deckname, Deck())


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
    """Useful decks for electrophysiology simulations."""

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
