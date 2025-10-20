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

# flake8: noqa: E501
import typing

from ansys.dyna.core.lib.card import Card, Field
from ansys.dyna.core.lib.keyword_base import KeywordBase


class EmEpCreatefiberorientation(KeywordBase):
    """LS-DYNA EM_EP_CREATEFIBERORIENTATION keyword"""

    keyword = "EM"
    subkeyword = "EP_CREATEFIBERORIENTATION"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field("partsid", int, 0, 10, kwargs.get("partsid")),
                    Field("solvid1", int, 10, 10, kwargs.get("solvid1")),
                    Field("solvid2", int, 20, 10, kwargs.get("solvid2", 1)),
                    Field("alpha", int, 30, 10, kwargs.get("alpha")),
                    Field("beta", int, 40, 10, kwargs.get("beta")),
                    Field("wfile", int, 50, 10, kwargs.get("wfile")),
                    Field("prerun", int, 60, 10, kwargs.get("prerun")),
                ],
            ),
        ]

    @property
    def partsid(self) -> typing.Optional[int]:
        """Get or set the Part set on which the system is solved"""  # nopep8
        return self._cards[0].get_value("partsid")

    @partsid.setter
    def partsid(self, value: int) -> None:
        self._cards[0].set_value("partsid", value)

    @property
    def solvid1(self) -> typing.Optional[int]:
        """Get or set the ID of the Laplace system that is solved in the transmural direction"""  # nopep8
        return self._cards[0].get_value("solvid1")

    @solvid1.setter
    def solvid1(self, value: int) -> None:
        self._cards[0].set_value("solvid1", value)

    @property
    def solvid2(self) -> int:
        """Get or set the ID of the Laplace system that is solved in the apicobasal direction"""  # nopep8
        return self._cards[0].get_value("solvid2")

    @solvid2.setter
    def solvid2(self, value: int) -> None:
        self._cards[0].set_value("solvid2", value)

    @property
    def alpha(self) -> typing.Optional[float]:
        """Get or set the helical angle with respect to the counterclockwise circumferential direction in the heart when looking from the base towards the apex. If a negative value is entered, a *DEFINE_â€ŒFUNCTION will be expected. See remark 1- for available parameters"""  # nopep8
        return self._cards[0].get_value("alpha")

    @alpha.setter
    def alpha(self, value: float) -> None:
        self._cards[0].set_value("alpha", value)

    @property
    def beta(self) -> typing.Optional[float]:
        """Get or set the angle with respect to the outward transmural axis of the heart. If a negative value is entered, a *DEFINE_â€ŒFUNCTION will be expected. See remark 1- for available parameters"""  # nopep8
        return self._cards[0].get_value("beta")

    @beta.setter
    def beta(self, value: float) -> None:
        self._cards[0].set_value("beta", value)

    @property
    def wfile(self) -> typing.Optional[int]:
        """Get or set the Selects whether result files (ELEMENT_SOLID_ORTHO.k and vtk files) are exported. Eq 0: not exported. Eq 1: exported"""  # nopep8
        return self._cards[0].get_value("wfile")

    @wfile.setter
    def wfile(self, value: int) -> None:
        self._cards[0].set_value("wfile", value)

    @property
    def prerun(self) -> typing.Optional[int]:
        """Get or set the Select whether the run is stopped after creating fibers: Eq. 0: do not stop after fiber creation. Eq 1: stop after fiber creation"""  # nopep8
        return self._cards[0].get_value("prerun")

    @prerun.setter
    def prerun(self, value: int) -> None:
        self._cards[0].set_value("prerun", value)
