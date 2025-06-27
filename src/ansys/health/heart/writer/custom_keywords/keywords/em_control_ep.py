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

from ansys.dyna.core.lib.card import Card, Field
from ansys.dyna.core.lib.keyword_base import KeywordBase


class EmControlEp(KeywordBase):
    """LS-DYNA EM_CONTROL_EP keyword"""

    keyword = "EM"
    subkeyword = "CONTROL_EP"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field("solvetype", int, 0, 10, kwargs.get("solvetype", 4)),
                    Field("numsplit", int, 10, 10, kwargs.get("numsplit", 1)),
                    Field("actusig", int, 20, 10, kwargs.get("actusig", 100000000)),
                    Field("ionsolvr", int, 30, 10, kwargs.get("ionsolvr")),
                ],
            ),
        ]

    @property
    def solvetype(self) -> int:
        """Get or set the ?"""  # nopep8
        return self._cards[0].get_value("solvetype")

    @solvetype.setter
    def solvetype(self, value: int) -> None:
        self._cards[0].set_value("solvetype", value)

    @property
    def numsplit(self) -> int:
        """Get or set the ?"""  # nopep8
        return self._cards[0].get_value("numsplit")

    @numsplit.setter
    def numsplit(self, value: int) -> None:
        self._cards[0].set_value("numsplit", value)

    @property
    def actusig(self) -> int:
        """Get or set the ?"""  # nopep8
        return self._cards[0].get_value("actusig")

    @actusig.setter
    def actusig(self, value: int) -> None:
        self._cards[0].set_value("actusig", value)

    @property
    def ionsolvr(self) -> int:
        """Get or set the ion solver type (0: Euler, 1: VODE, 2: Spline"""  # nopep8
        return self._cards[0].get_value("ionsolvr")

    @ionsolvr.setter
    def ionsolvr(self, value: int) -> None:
        self._cards[0].set_value("ionsolvr", value)
