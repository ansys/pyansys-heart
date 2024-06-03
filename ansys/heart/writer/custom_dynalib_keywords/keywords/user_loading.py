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

from ansys.dyna.keywords.lib.card import Card, Field
from ansys.dyna.keywords.lib.keyword_base import KeywordBase


class UserLoading(KeywordBase):
    """DYNA USER_LOADING keyword"""

    keyword = "USER"
    subkeyword = "LOADING"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field("parm1", float, 0, 10, kwargs.get("parm1")),
                    Field("parm2", float, 10, 10, kwargs.get("parm2")),
                    Field("parm3", float, 20, 10, kwargs.get("parm3")),
                    Field("parm4", float, 30, 10, kwargs.get("parm4")),
                    Field("parm5", float, 40, 10, kwargs.get("parm5")),
                    Field("parm6", float, 50, 10, kwargs.get("parm6")),
                    Field("parm7", float, 60, 10, kwargs.get("parm7")),
                    Field("parm8", float, 70, 10, kwargs.get("parm8")),
                ],
            ),
        ]

    @property
    def parm1(self) -> typing.Optional[float]:
        """Get or set the This is the Nth user input parameter"""  # nopep8
        return self._cards[0].get_value("parm1")

    @parm1.setter
    def parm1(self, value: float) -> None:
        self._cards[0].set_value("parm1", value)

    @property
    def parm2(self) -> typing.Optional[float]:
        """Get or set the This is the Nth user input parameter"""  # nopep8
        return self._cards[0].get_value("parm2")

    @parm2.setter
    def parm2(self, value: float) -> None:
        self._cards[0].set_value("parm2", value)

    @property
    def parm3(self) -> typing.Optional[float]:
        """Get or set the This is the Nth user input parameter"""  # nopep8
        return self._cards[0].get_value("parm3")

    @parm3.setter
    def parm3(self, value: float) -> None:
        self._cards[0].set_value("parm3", value)

    @property
    def parm4(self) -> typing.Optional[float]:
        """Get or set the This is the Nth user input parameter"""  # nopep8
        return self._cards[0].get_value("parm4")

    @parm4.setter
    def parm4(self, value: float) -> None:
        self._cards[0].set_value("parm4", value)

    @property
    def parm5(self) -> typing.Optional[float]:
        """Get or set the This is the Nth user input parameter"""  # nopep8
        return self._cards[0].get_value("parm5")

    @parm5.setter
    def parm5(self, value: float) -> None:
        self._cards[0].set_value("parm5", value)

    @property
    def parm6(self) -> typing.Optional[float]:
        """Get or set the This is the Nth user input parameter"""  # nopep8
        return self._cards[0].get_value("parm6")

    @parm6.setter
    def parm6(self, value: float) -> None:
        self._cards[0].set_value("parm6", value)

    @property
    def parm7(self) -> typing.Optional[float]:
        """Get or set the This is the Nth user input parameter"""  # nopep8
        return self._cards[0].get_value("parm7")

    @parm7.setter
    def parm7(self, value: float) -> None:
        self._cards[0].set_value("parm7", value)

    @property
    def parm8(self) -> typing.Optional[float]:
        """Get or set the This is the Nth user input parameter"""  # nopep8
        return self._cards[0].get_value("parm8")

    @parm8.setter
    def parm8(self, value: float) -> None:
        self._cards[0].set_value("parm8", value)
