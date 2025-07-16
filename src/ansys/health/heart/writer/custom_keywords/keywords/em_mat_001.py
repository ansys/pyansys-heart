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


class EmMat001(KeywordBase):
    """LS-DYNA EM_MAT_001 keyword"""

    keyword = "EM"
    subkeyword = "MAT_001"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field("mid", int, 0, 10, kwargs.get("mid")),
                    Field("mtype", int, 10, 10, kwargs.get("mtype", 0)),
                    Field("sigma", float, 20, 10, kwargs.get("sigma")),
                    Field("eosid", int, 30, 10, kwargs.get("eosid")),
                    Field("unused", int, 40, 10, kwargs.get("unused")),
                    Field("lambda_", float, 50, 10, kwargs.get("lambda_")),
                    Field("deatht", float, 60, 10, kwargs.get("deatht", 1e28)),
                    Field("rdltype", int, 70, 10, kwargs.get("rdltype", 0)),
                ],
            ),
            Card(
                [
                    Field("empty", float, 0, 10, kwargs.get("empty", None)),
                    Field("freq", float, 10, 10, kwargs.get("freq", None)),
                    Field("epsrr", float, 20, 10, kwargs.get("epsrr", None)),
                    Field("eosid2", int, 30, 10, kwargs.get("eosid2", None)),
                    Field("epsri", float, 40, 10, kwargs.get("epsri", None)),
                    Field("eosid3", int, 50, 10, kwargs.get("eosid3", None)),
                    Field("eps0", float, 60, 10, kwargs.get("eps0", None)),
                ],
            ),
            Card(
                [
                    Field("beta", float, 0, 10, kwargs.get("beta", 140)),
                    Field("cm", float, 10, 10, kwargs.get("cm", 0.01)),
                ],
            ),
        ]

    @property
    def mid(self) -> typing.Optional[int]:
        """Get or set the Material ID: refers to MID in the *PART card."""  # nopep8
        return self._cards[0].get_value("mid")

    @mid.setter
    def mid(self, value: int) -> None:
        self._cards[0].set_value("mid", value)

    @property
    def mtype(self) -> int:
        """Get or set the defines the electromagnetism type of the material:
        EQ.0: Air or vacuum
        EQ.1: Insulator material: these materials have the same electromagnetism behavior as EQ.0
        EQ.2: Conductor carrying a source. In these conductors, the eddy current problem is solved,
            which gives the actual current density. Typically, this would correspond to the coil.
        EQ.3: Fluid conductor. In that case, MID refers to the ID given in * ICFD_PART.
        EQ.4: Conductor not connected to any current or voltage source, where the Eddy
            current problem is solved. Typically, this would correspond to the workpiece.
        """  # nopep8
        return self._cards[0].get_value("mtype")

    @mtype.setter
    def mtype(self, value: int) -> None:
        if value not in [0, 1, 2, 3, 4]:
            raise Exception("""mtype must be one of {0,1,2,3,4}""")
        self._cards[0].set_value("mtype", value)

    @property
    def sigma(self) -> typing.Optional[float]:
        """Get or set the initial electrical conductivity of the material."""  # nopep8
        return self._cards[0].get_value("sigma")

    @sigma.setter
    def sigma(self, value: float) -> None:
        self._cards[0].set_value("sigma", value)

    @property
    def eosid(self) -> typing.Optional[int]:
        """Get or set the Optional ID of the EOS to be used for the electrical conductivity
        (see *EM_EOS card)."""  # nopep8
        return self._cards[0].get_value("eosid")

    @eosid.setter
    def eosid(self, value: int) -> None:
        self._cards[0].set_value("eosid", value)

    @property
    def lambda_(self) -> float:
        """Get or set the conductivity ratio value lambda."""  # nopep8
        return self._cards[0].get_value("deatht")

    @lambda_.setter
    def lambda_(self, value: float) -> None:
        self._cards[0].set_value("lambda_", value)

    @property
    def deatht(self) -> float:
        """Get or set the Death time for the material. After DEATHT, the material will no longer be
        considered a conductor and removed from the EM solve."""  # nopep8
        return self._cards[0].get_value("deatht")

    @deatht.setter
    def deatht(self, value: float) -> None:
        self._cards[0].set_value("deatht", value)

    @property
    def rdltype(self) -> int:
        """Get or set the Used for the application: composite Tshell battery, with
        **EM_RANDLES_TSHELL.Defines the function of the layer associated to MID:
        EQ.0:	Default. Conductor which is not part of a battery cell
        EQ.1:Current Collector Positive
        EQ.2: Positive Electrode
        EQ.3:Separator
        EQ.4:Negative Electrode
        EQ.5:Current Collector Negative
        """  # nopep8
        return self._cards[0].get_value("rdltype")

    @rdltype.setter
    def rdltype(self, value: int) -> None:
        if value not in [0, 1, 2, 3, 4, 5]:
            raise Exception("""rdltype must be one of {0,1,2,3,4,5}""")
        self._cards[0].set_value("rdltype", value)

    @property
    def freq(self) -> typing.Optional[float]:
        """freq."""  # nopep8
        return self._cards[1].get_value("freq")

    @freq.setter
    def freq(self, value: float) -> None:
        self._cards[1].set_value("freq", value)

    @property
    def epsrr(self) -> typing.Optional[float]:
        """epsri."""  # nopep8
        return self._cards[1].get_value("epsrr")

    @epsrr.setter
    def epsrr(self, value: float) -> None:
        self._cards[1].set_value("epsrr", value)

    @property
    def eosid2(self) -> typing.Optional[int]:
        """eosid2."""  # nopep8
        return self._cards[1].get_value("eosid2")

    @eosid2.setter
    def eosid2(self, value: int) -> None:
        self._cards[1].set_value("eosid2", value)

    @property
    def epsri(self) -> typing.Optional[float]:
        """epsri."""  # nopep8
        return self._cards[1].get_value("epsri")

    @epsri.setter
    def epsri(self, value: float) -> None:
        self._cards[1].set_value("epsri", value)

    @property
    def eosid3(self) -> typing.Optional[int]:
        """eosid3."""  # nopep8
        return self._cards[1].get_value("eosid3")

    @eosid3.setter
    def eosid3(self, value: int) -> None:
        self._cards[1].set_value("eosid3", value)

    @property
    def eps0(self) -> typing.Optional[float]:
        """eps0."""  # nopep8
        return self._cards[1].get_value("eps0")

    @eps0.setter
    def eps0(self, value: float) -> None:
        self._cards[1].set_value("eps0", value)

    @property
    def beta(self) -> typing.Optional[float]:
        """Surface to volume ratio."""  # nopep8
        return self._cards[2].get_value("beta")

    @beta.setter
    def beta(self, value: float) -> None:
        self._cards[2].set_value("beta", value)

    @property
    def cm(self) -> typing.Optional[float]:
        """Membrane capacitance."""  # nopep8
        return self._cards[2].get_value("cm")

    @cm.setter
    def cm(self, value: float) -> None:
        self._cards[2].set_value("cm", value)
