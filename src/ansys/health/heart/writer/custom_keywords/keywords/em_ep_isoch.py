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

import typing

from ansys.dyna.core.lib.card import Card, Field
from ansys.dyna.core.lib.keyword_base import KeywordBase


class EmEpIsoch(KeywordBase):
    """LS-DYNA EM_EP_ISOCH keyword"""

    keyword = "EM"
    subkeyword = "EP_ISOCH"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field("idisoch", int, 0, 10, kwargs.get("idisoch")),
                    Field("idepol", int, 10, 10, kwargs.get("idepol", 0)),
                    Field("dplthr", float, 20, 10, kwargs.get("dplthr", 0)),
                    Field("irepol", int, 30, 10, kwargs.get("irepol", 0)),
                    Field("rplthr", float, 40, 10, kwargs.get("rplthr", 0)),
                    Field("cyclenmin", float, 50, 10, kwargs.get("cyclenmin", 2000)),
                    Field("apdmin", float, 60, 10, kwargs.get("apdmin", 100)),
                ],
            ),
        ]

    @property
    def idisoch(self) -> typing.Optional[int]:
        """Get or set the ID of the isochrone."""  # nopep8
        return self._cards[0].get_value("idisoch")

    @idisoch.setter
    def idisoch(self, value: int) -> None:
        self._cards[0].set_value("idisoch", value)

    @property
    def idepol(self) -> int:
        """Get or set the Flag to activate the computation of depolarization:
        EQ.0: OFF
        EQ.1:ON
        """  # nopep8
        return self._cards[0].get_value("idepol")

    @idepol.setter
    def idepol(self, value: int) -> None:
        self._cards[0].set_value("idepol", value)

    @property
    def dplthr(self) -> float:
        """Get or set the Amplitude threshold used for measuring depolarization."""  # nopep8
        return self._cards[0].get_value("dplthr")

    @dplthr.setter
    def dplthr(self, value: float) -> None:
        self._cards[0].set_value("dplthr", value)

    @property
    def irepol(self) -> int:
        """Get or set the Flag to activate the computation of repolarization:
        EQ.0: OFF
        EQ.1:ON
        """  # nopep8
        return self._cards[0].get_value("irepol")

    @irepol.setter
    def irepol(self, value: int) -> None:
        self._cards[0].set_value("irepol", value)

    @property
    def rplthr(self) -> float:
        """Get or set the Amplitude threshold used for measuring repolarization."""  # nopep8
        return self._cards[0].get_value("rplthr")

    @rplthr.setter
    def rplthr(self, value: float) -> None:
        self._cards[0].set_value("rplthr", value)

    @property
    def cyclenmin(self) -> float:
        """Get or set the minimum delay between two repol and two depol."""  # nopep8
        return self._cards[0].get_value("cyclenmin")

    @cyclenmin.setter
    def cyclenmin(self, value: float) -> None:
        self._cards[0].set_value("cyclenmin", value)

    @property
    def apdmin(self) -> float:
        """Get or set the minimum value to be considered for APD."""  # nopep8
        return self._cards[0].get_value("apdmin")

    @apdmin.setter
    def apdmin(self, value: float) -> None:
        self._cards[0].set_value("apdmin", value)
