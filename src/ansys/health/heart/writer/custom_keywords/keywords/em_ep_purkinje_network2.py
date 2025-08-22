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


class EmEpPurkinjeNetwork2(KeywordBase):
    """LS-DYNA EM_EP_PURKINJE_NETWORK keyword"""

    keyword = "EM"
    subkeyword = "EP_PURKINJE_NETWORK2"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field("purkid", int, 0, 10, kwargs.get("purkid")),
                    Field("buildnet", int, 10, 10, kwargs.get("buildnet")),
                    Field("ssid", int, 20, 10, kwargs.get("ssid")),
                    Field("mid", int, 30, 10, kwargs.get("mid")),
                    Field("nodestart", float, 40, 10, kwargs.get("nodestart")),
                    Field("pointstx", float, 50, 10, kwargs.get("pointstx")),
                    Field("pointsty", float, 60, 10, kwargs.get("pointsty")),
                    Field("pointstz", float, 70, 10, kwargs.get("pointstz")),
                ],
            ),
            Card(
                [
                    Field("edgelen", float, 0, 10, kwargs.get("edgelen", 3)),
                    Field("ngen", int, 10, 10, kwargs.get("ngen", 200)),
                    Field("nbrinit", int, 20, 10, kwargs.get("nbrinit", 8)),
                    Field("nsplit", int, 30, 10, kwargs.get("nsplit", 2)),
                    Field("inodeid", int, 40, 10, kwargs.get("inodeid")),
                    Field("iedgeid", int, 50, 10, kwargs.get("iedgeid")),
                ],
            ),
            Card(
                [
                    Field("pmjtype", int, 0, 10, kwargs.get("pmjtype", 1)),
                    Field("pmjradius", float, 10, 10, kwargs.get("pmjradius", 0.7)),
                    Field("pmjrestype", int, 20, 10, kwargs.get("pmjrestype", 1)),
                    Field("pmjres", float, 30, 10, kwargs.get("pmjres", 1)),
                ],
            ),
        ]

    @property
    def purkid(self) -> typing.Optional[int]:
        """Get or set the ID for the Purkinje network"""  # nopep8
        return self._cards[0].get_value("purkid")

    @purkid.setter
    def purkid(self, value: int) -> None:
        self._cards[0].set_value("purkid", value)

    @property
    def buildnet(self) -> typing.Optional[int]:
        """Get or set the Flag to create Purkinje network: EQ.0: Purkinje network not created, EQ.1: New Purkinje network created."""  # nopep8
        return self._cards[0].get_value("buildnet")

    @buildnet.setter
    def buildnet(self, value: int) -> None:
        self._cards[0].set_value("buildnet", value)

    @property
    def ssid(self) -> typing.Optional[int]:
        """Get or set the Segment set on which the Purkinje network is lying"""  # nopep8
        return self._cards[0].get_value("ssid")

    @ssid.setter
    def ssid(self, value: int) -> None:
        self._cards[0].set_value("ssid", value)

    @property
    def mid(self) -> typing.Optional[int]:
        """Get or set the Material ID defined in the *MAT section"""  # nopep8
        return self._cards[0].get_value("mid")

    @mid.setter
    def mid(self, value: int) -> None:
        self._cards[0].set_value("mid", value)

    @property
    def nodestart(self) -> typing.Optional[int]:
        """Get or set the Material ID defined in the *MAT section"""  # nopep8
        return self._cards[0].get_value("nodestart")

    @nodestart.setter
    def nodestart(self, value: int) -> None:
        self._cards[0].set_value("nodestart", value)

    @property
    def pointstx(self) -> typing.Optional[float]:
        """Get or set the X coordinate of the tree origin"""  # nopep8
        return self._cards[0].get_value("pointstx")

    @pointstx.setter
    def pointstx(self, value: float) -> None:
        self._cards[0].set_value("pointstx", value)

    @property
    def pointsty(self) -> typing.Optional[float]:
        """Get or set the Y coordinate of the tree origin"""  # nopep8
        return self._cards[0].get_value("pointsty")

    @pointsty.setter
    def pointsty(self, value: float) -> None:
        self._cards[0].set_value("pointsty", value)

    @property
    def pointstz(self) -> typing.Optional[float]:
        """Get or set the Z coordinate of the tree origin"""  # nopep8
        return self._cards[0].get_value("pointstz")

    @pointstz.setter
    def pointstz(self, value: float) -> None:
        self._cards[0].set_value("pointstz", value)

    @property
    def edgelen(self) -> float:
        """Get or set the Edge length"""  # nopep8
        return self._cards[1].get_value("edgelen")

    @edgelen.setter
    def edgelen(self, value: float) -> None:
        self._cards[1].set_value("edgelen", value)

    @property
    def ngen(self) -> int:
        """Get or set the Number of generations of branches"""  # nopep8
        return self._cards[1].get_value("ngen")

    @ngen.setter
    def ngen(self, value: int) -> None:
        self._cards[1].set_value("ngen", value)

    @property
    def nbrinit(self) -> int:
        """Get or set the Number of branches attached to the tree origin"""  # nopep8
        return self._cards[1].get_value("nbrinit")

    @nbrinit.setter
    def nbrinit(self, value: int) -> None:
        self._cards[1].set_value("nbrinit", value)

    @property
    def nsplit(self) -> int:
        """Get or set the Number of child branches at each node of the tree"""  # nopep8
        return self._cards[1].get_value("nsplit")

    @nsplit.setter
    def nsplit(self, value: int) -> None:
        self._cards[1].set_value("nsplit", value)

    @property
    def inodeid(self) -> typing.Optional[int]:
        """Get or set the Initial node ID"""  # nopep8
        return self._cards[1].get_value("inodeid")

    @inodeid.setter
    def inodeid(self, value: int) -> None:
        self._cards[1].set_value("inodeid", value)

    @property
    def iedgeid(self) -> typing.Optional[int]:
        """Get or set the Initial edge ID"""  # nopep8
        return self._cards[1].get_value("iedgeid")

    @iedgeid.setter
    def iedgeid(self, value: int) -> None:
        self._cards[1].set_value("iedgeid", value)

    @property
    def pmjtype(self) -> typing.Optional[int]:
        """Purkinje muscle junction type:
        Eq. 0: The resistance of the beam is given by R = 1/sigma * L where L is the beam length.
        Eq. 1: Connect beam leaves to the closest solid node.
        Eq. 2: Connect beam leaves with solid nodes within a disk.
        Eq. 3: Connect all beam nodes with the closest solid nodes.
        Eq. 4: Connect all beam nodes with the solid nodes within a disk.
        """  # nopep8
        return self._cards[2].get_value("pmjtype")

    @pmjtype.setter
    def pmjtype(self, value: int) -> None:
        self._cards[2].set_value("pmjtype", value)

    @property
    def pmjradius(self) -> typing.Optional[float]:
        """Radius of the disk connecting beam nodes with solid nodes."""  # nopep8
        return self._cards[2].get_value("pmjradius")

    @pmjradius.setter
    def pmjradius(self, value: float) -> None:
        self._cards[2].set_value("pmjradius", value)

    @property
    def pmjrestype(self) -> typing.Optional[int]:
        """Purkinje muscle junction resistance type:
        Eq. 1: Apply resistance to each beam.
        """  # nopep8
        return self._cards[2].get_value("pmjrestype")

    @pmjrestype.setter
    def pmjrestype(self, value: int) -> None:
        self._cards[2].set_value("pmjrestype", value)

    @property
    def pmjres(self) -> typing.Optional[float]:
        """Beam resistance value."""  # nopep8
        return self._cards[2].get_value("pmjres")

    @pmjres.setter
    def pmjres(self, value: float) -> None:
        self._cards[2].set_value("pmjres", value)
