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


class EmEpCellmodelTomek(KeywordBase):
    """LS-DYNA EM_EP_CELLMODEL_TOMEK keyword"""

    keyword = "EM"
    subkeyword = "EP_CELLMODEL_TOMEK"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field("mid", int, 0, 10, kwargs.get("mid", 1)),
                    Field("phiendmid", float, 10, 10, kwargs.get("phiendmid", 0.17)),
                    Field("phimidepi", float, 20, 10, kwargs.get("phimidepi", 0.58)),
                ],
            ),
        ]

    @property
    def mid(self) -> int:
        """Get or set the Material ID"""  # nopep8
        return self._cards[0].get_value("mid")

    @mid.setter
    def mid(self, value: int) -> None:
        self._cards[0].set_value("mid", value)

    @property
    def phiendmid(self) -> float:
        """Get or set the Phi endocardium > mid"""  # nopep8
        return self._cards[0].get_value("phiendmid")

    @phiendmid.setter
    def phiendmid(self, value: float) -> None:
        self._cards[0].set_value("phiendmid", value)

    @property
    def phimidepi(self) -> float:
        """Get or set the Phi mid > epicardium"""  # nopep8
        return self._cards[0].get_value("phimidepi")

    @phimidepi.setter
    def phimidepi(self, value: float) -> None:
        self._cards[0].set_value("phimidepi", value)
