# Copyright (C) 2023 - 2024 ANSYS, Inc. and/or its affiliates.
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
from ansys.dyna.core.lib.duplicate_card import DuplicateCard
from ansys.dyna.core.lib.keyword_base import KeywordBase


class ConstrainedLinearGlobal(KeywordBase):
    """DYNA CONSTRAINED_LINEAR_GLOBAL keyword"""

    keyword = "CONSTRAINED"
    subkeyword = "LINEAR_GLOBAL"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field("licd", int, 0, 10, kwargs.get("licd")),
                ],
            ),
            DuplicateCard(
                [
                    Field("nid", int, 0, 10, kwargs.get("nid")),
                    Field("dof", int, 10, 10, kwargs.get("dof")),
                    Field("coef", float, 20, 10, kwargs.get("coef")),
                ],
                None,
            ),
        ]

    @property
    def linear_constraints(self):
        """Gets the table of linear constraints"""
        return self._cards[1].table

    @linear_constraints.setter
    def linear_constraints(self, df):
        """sets linear constraints from the dataframe df"""
        self._cards[1].table = df

    @property
    def licd(self) -> typing.Optional[int]:
        """Get or set the Linear constraint definition ID. This ID can be used
        to identify a set to which this constraint is a member."""  # nopep8
        return self._cards[0].get_value("licd")

    @licd.setter
    def licd(self, value: int) -> None:
        self._cards[0].set_value("licd", value)
