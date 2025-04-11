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


class EmEpTentusscherStimulus2(KeywordBase):
    """DYNA EM_EP_TENTUSSCHER_STIMULUS2 keyword"""

    keyword = "EM"
    subkeyword = "EP_TENTUSSCHER_STIMULUS2"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field("stimid", int, 0, 10, kwargs.get("stimid")),
                    Field("settype", int, 10, 10, kwargs.get("settype")),
                    Field("setid", int, 20, 10, kwargs.get("setid")),
                    Field("lcid", int, 30, 10, kwargs.get("lcid")),
                ],
            ),
        ]

    @property
    def stimid(self) -> typing.Optional[int]:
        """Get or set the ID of the stimulation"""  # nopep8
        return self._cards[0].get_value("stimid")

    @stimid.setter
    def stimid(self, value: int) -> None:
        self._cards[0].set_value("stimid", value)

    @property
    def settype(self) -> typing.Optional[int]:
        """Get or set the Set type: EQ.1: Segment set, EQ.2: Node set"""  # nopep8
        return self._cards[0].get_value("settype")

    @settype.setter
    def settype(self, value: int) -> None:
        self._cards[0].set_value("settype", value)

    @property
    def setid(self) -> typing.Optional[int]:
        """Get or set the Node set or segment set ID to be stimulated"""  # nopep8
        return self._cards[0].get_value("setid")

    @setid.setter
    def setid(self, value: int) -> None:
        self._cards[0].set_value("setid", value)
    @property
    def lcid(self) -> typing.Optional[int]:
        """Load curve id defined to create the stimulation"""  # nopep8
        return self._cards[0].get_value("lcid")

    @lcid.setter
    def lcid(self, value: int) -> None:
        self._cards[0].set_value("lcid", value)
