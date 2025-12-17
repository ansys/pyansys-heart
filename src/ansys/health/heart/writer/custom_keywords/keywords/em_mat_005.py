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


class EmMat005(KeywordBase):
    """LS-DYNA EM_MAT_005 keyword"""

    keyword = "EM"
    subkeyword = "MAT_005"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field("mid", int, 0, 10, kwargs.get("mid")),
                    Field("mtype", int, 10, 10, kwargs.get("mtype", 0)),
                    Field("cvsigma11", float, 20, 10, kwargs.get("cvsigma11")),
                    Field("cvsigma22", float, 30, 10, kwargs.get("cvsigma22")),
                    Field("cvsigma33", float, 40, 10, kwargs.get("cvsigma33")),
                    Field("beta", float, 50, 10, kwargs.get("beta", 0.14)),
                    Field("cm", float, 60, 10, kwargs.get("cm", 1.0)),
                ],
            ),
            Card(
                [
                    Field("sigma12", float, 0, 10, kwargs.get("sigma12")),
                    Field("sigma13", float, 10, 10, kwargs.get("sigma13")),
                    Field("sigma21", float, 20, 10, kwargs.get("sigma21")),
                    Field("sigma23", float, 30, 10, kwargs.get("sigma23")),
                    Field("sigma31", float, 40, 10, kwargs.get("sigma31")),
                    Field("sigma32", float, 50, 10, kwargs.get("sigma32")),
                ],
            ),
            Card(
                [
                    Field("condsigma11", float, 20, 10, kwargs.get("condsigma11")),
                    Field("condsigma22", float, 30, 10, kwargs.get("condsigma22")),
                    Field("condsigma33", float, 40, 10, kwargs.get("condsigma33")),
                ],
            ),
            Card(
                [
                    Field("AOPT", int, 0, 10, kwargs.get("AOPT",2)),
                    Field("xp", float, 10, 10, kwargs.get("xp")),
                    Field("yp", float, 20, 10, kwargs.get("yp")),
                    Field("zp", float, 30, 10, kwargs.get("zp")),
                    Field("a1", float, 40, 10, kwargs.get("a1")),
                    Field("a2", float, 50, 10, kwargs.get("a2")),
                    Field("a3", float, 60, 10, kwargs.get("a3")),
                    Field("macf", float, 70, 10, kwargs.get("macf",1)),
                ],
            ),
            Card(
                [
                    Field("v1", float, 0, 10, kwargs.get("v1")),
                    Field("v2", float, 10, 10, kwargs.get("v2")),
                    Field("v3", float, 20, 10, kwargs.get("v3")),
                    Field("d1", float, 30, 10, kwargs.get("d1")),
                    Field("d2", float, 40, 10, kwargs.get("d2")),
                    Field("d3", float, 50, 10, kwargs.get("d3")),
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
        """Get or set the Defines the electromagnetism type of the material:
            EQ.0: Air or vacuum
            EQ.1: Insulator material: these materials have the same electromagnetism behavior as EQ.0
            EQ.2: Conductor carrying a source. In these conductors, the eddy current problem is solved, which gives the actual current density. Typically, this would correspond to the coil.
            EQ.4: Conductor not connected to any current or voltage source, where the Eddy current problem is solved. Typically, this would correspond to the workpiece
        .
        """  # nopep8
        return self._cards[0].get_value("mtype")

    @mtype.setter
    def mtype(self, value: int) -> None:
        if value not in [0, 1, 2, 4]:
            raise Exception("""mtype must be one of {0,1,2,4}""")
        self._cards[0].set_value("mtype", value)

    @property
    def CVsigma11(self) -> typing.Optional[float]:
        """Get or set the The 1,1 term in the 3 x 3 electromagnetic conductivity tensor matrix. Note that 1 corresponds to the a material direction.If a negative value is entered, a *DEFINE_FUNCTION will be expected. See remark 3- for available parameters."""  # nopep8
        return self._cards[0].get_value("CVsigma11")

    @CVsigma11.setter
    def CVsigma11(self, value: float) -> None:
        self._cards[0].set_value("CVsigma11", value)

    @property
    def CVsigma22(self) -> typing.Optional[float]:
        """Get or set the The 2,2 term in the 3 x 3 electromagnetic conductivity tensor matrix.If a negative value is entered, a *DEFINE_FUNCTION will be expected. See remark 3- for available parameters."""  # nopep8
        return self._cards[0].get_value("CVsigma22")

    @CVsigma22.setter
    def CVsigma22(self, value: float) -> None:
        self._cards[0].set_value("CVsigma22", value)

    @property
    def CVsigma33(self) -> typing.Optional[float]:
        """Get or set the The 2,2 term in the 3 x 3 electromagnetic conductivity tensor matrix.If a negative value is entered, a *DEFINE_FUNCTION will be expected. See remark 3- for available parameters."""  # nopep8
        return self._cards[0].get_value("CVsigma33")

    @CVsigma33.setter
    def CVsigma33(self, value: float) -> None:
        self._cards[0].set_value("CVsigma33", value)

    @property
    def beta(self) -> float:
        """Get or set the Surface to volume ratio."""  # nopep8
        return self._cards[0].get_value("beta")

    @beta.setter
    def beta(self, value: float) -> None:
        self._cards[0].set_value("beta", value)

    @property
    def cm(self) -> float:
        """Get or set the Membrane capacitance"""  # nopep8
        return self._cards[0].get_value("cm")

    @cm.setter
    def cm(self, value: float) -> None:
        self._cards[0].set_value("cm", value)

    @property
    def sigma12(self) -> typing.Optional[float]:
        """Get or set the The 1,2 term in the 3 x 3 electromagnetic conductivity tensor matrix.Note that 2 corresponds to the b material direction.. If a negative value is entered, a *DEFINE_FUNCTION will be expected. See remark 3- for available parameters."""  # nopep8
        return self._cards[1].get_value("sigma12")

    @sigma12.setter
    def sigma12(self, value: float) -> None:
        self._cards[1].set_value("sigma12", value)

    @property
    def sigma13(self) -> typing.Optional[float]:
        """Get or set the The 1,3 term in the 3 x 3 electromagnetic conductivity tensor matrix.If a negative value is entered, a *DEFINE_FUNCTION will be expected. See remark 3- for available parameters."""  # nopep8
        return self._cards[1].get_value("sigma13")

    @sigma13.setter
    def sigma13(self, value: float) -> None:
        self._cards[1].set_value("sigma13", value)

    @property
    def sigma21(self) -> typing.Optional[float]:
        """Get or set the The 2,1 term in the 3 x 3 electromagnetic conductivity tensor matrix. Note that 1 corresponds to the a material direction.If a negative value is entered, a *DEFINE_FUNCTION will be expected. See remark 3- for available parameters."""  # nopep8
        return self._cards[1].get_value("sigma21")

    @sigma21.setter
    def sigma21(self, value: float) -> None:
        self._cards[1].set_value("sigma21", value)

    @property
    def sigma23(self) -> typing.Optional[float]:
        """Get or set the The 2,3 term in the 3 x 3 electromagnetic conductivity tensor matrix.If a negative value is entered, a *DEFINE_FUNCTION will be expected. See remark 3- for available parameters."""  # nopep8
        return self._cards[1].get_value("sigma23")

    @sigma23.setter
    def sigma23(self, value: float) -> None:
        self._cards[1].set_value("sigma23", value)

    @property
    def sigma31(self) -> typing.Optional[float]:
        """Get or set the The 3,1 term in the 3 x 3 electromagnetic conductivity tensor matrix.If a negative value is entered, a *DEFINE_FUNCTION will be expected. See remark 3- for available parameters."""  # nopep8
        return self._cards[1].get_value("sigma31")

    @sigma31.setter
    def sigma31(self, value: float) -> None:
        self._cards[1].set_value("sigma31", value)

    @property
    def sigma32(self) -> typing.Optional[float]:
        """Get or set the The 3,2 term in the 3 x 3 electromagnetic conductivity tensor matrix.If a negative value is entered, a *DEFINE_FUNCTION will be expected. See remark 3- for available parameters."""  # nopep8
        return self._cards[1].get_value("sigma32")

    @sigma32.setter
    def sigma32(self, value: float) -> None:
        self._cards[1].set_value("sigma32", value)


    @property
    def Condsigma11(self) -> typing.Optional[float]:
          return self._cards[2].get_value("Condsigma11")

    @sigma32.setter
    def Condsigma11(self, value: float) -> None:
        self._cards[2].set_value("Condsigma11", value)

    @property
    def Condsigma22(self) -> typing.Optional[float]:
          return self._cards[2].get_value("Condsigma22")

    @sigma32.setter
    def Condsigma22(self, value: float) -> None:
        self._cards[2].set_value("Condsigma22", value)

    @property
    def Condsigma33(self) -> typing.Optional[float]:
          return self._cards[2].get_value("Condsigma33")

    @sigma32.setter
    def Condsigma33(self, value: float) -> None:
        self._cards[2].set_value("Condsigma33", value)

    @property
    def xp1(self) -> typing.Optional[float]:
        """Get or set the Define coordinates of point p for AOPT = 1 and 4."""  # nopep8
        return self._cards[2].get_value("xp1")

    @xp1.setter
    def xp1(self, value: float) -> None:
        self._cards[2].set_value("xp1", value)

    @property
    def yp1(self) -> typing.Optional[float]:
        """Get or set the Define coordinates of point p for AOPT = 1 and 4."""  # nopep8
        return self._cards[2].get_value("yp1")

    @yp1.setter
    def yp1(self, value: float) -> None:
        self._cards[2].set_value("yp1", value)

    @property
    def xp(self) -> typing.Optional[float]:
        """Get or set the Define coordinates of point p for AOPT = 1 and 4."""  # nopep8
        return self._cards[3].get_value("xp")

    @xp.setter
    def xp(self, value: float) -> None:
        self._cards[3].set_value("xp", value)

    @property
    def yp(self) -> typing.Optional[float]:
        """Get or set the Define coordinates of point p for AOPT = 1 and 4."""  # nopep8
        return self._cards[3].get_value("yp")

    @yp.setter
    def yp(self, value: float) -> None:
        self._cards[3].set_value("yp", value)

    @property
    def zp(self) -> typing.Optional[float]:
        """Get or set the Define coordinates of point p for AOPT = 1 and 4."""  # nopep8
        return self._cards[3].get_value("zp")

    @zp.setter
    def zp(self, value: float) -> None:
        self._cards[3].set_value("zp", value)

    @property
    def a1(self) -> typing.Optional[float]:
        """Get or set the Define components of vector a for AOPT = 2."""  # nopep8
        return self._cards[3].get_value("a1")

    @a1.setter
    def a1(self, value: float) -> None:
        self._cards[3].set_value("a1", value)

    @property
    def a2(self) -> typing.Optional[float]:
        """Get or set the Define components of vector a for AOPT = 2."""  # nopep8
        return self._cards[3].get_value("a2")

    @a2.setter
    def a2(self, value: float) -> None:
        self._cards[3].set_value("a2", value)

    @property
    def a3(self) -> typing.Optional[float]:
        """Get or set the Define components of vector a for AOPT = 2."""  # nopep8
        return self._cards[3].get_value("a3")

    @a3.setter
    def a3(self, value: float) -> None:
        self._cards[3].set_value("a3", value)

    @property
    def macf(self) -> int:
        """Get or set the Material axes change flag for solid elements:
        EQ.1: No change, default
        """  # nopep8
        return self._cards[3].get_value("macf")

    @macf.setter
    def macf(self, value: int) -> None:
        self._cards[3].set_value("macf", value)

    @property
    def AOPT(self) -> int:
        return self._cards[3].get_value("AOPT")

    @AOPT.setter
    def AOPT(self, value: int) -> None:
        self._cards[3].set_value("AOPT", value)


    @property
    def v1(self) -> typing.Optional[float]:
        """Get or set the Define components of vector v for AOPT = 3 and 4."""  # nopep8
        return self._cards[4].get_value("v1")

    @v1.setter
    def v1(self, value: float) -> None:
        self._cards[4].set_value("v1", value)

    @property
    def v2(self) -> typing.Optional[float]:
        """Get or set the Define components of vector v for AOPT = 3 and 4."""  # nopep8
        return self._cards[4].get_value("v2")

    @v2.setter
    def v2(self, value: float) -> None:
        self._cards[4].set_value("v2", value)

    @property
    def v3(self) -> typing.Optional[float]:
        """Get or set the Define components of vector v for AOPT = 3 and 4."""  # nopep8
        return self._cards[4].get_value("v3")

    @v3.setter
    def v3(self, value: float) -> None:
        self._cards[4].set_value("v3", value)

    @property
    def d1(self) -> typing.Optional[float]:
        """Get or set the Define components of vector d for AOPT = 2."""  # nopep8
        return self._cards[4].get_value("d1")

    @d1.setter
    def d1(self, value: float) -> None:
        self._cards[4].set_value("d1", value)

    @property
    def d2(self) -> typing.Optional[float]:
        """Get or set the Define components of vector d for AOPT = 2."""  # nopep8
        return self._cards[4].get_value("d2")

    @d2.setter
    def d2(self, value: float) -> None:
        self._cards[4].set_value("d2", value)

    @property
    def d3(self) -> typing.Optional[float]:
        """Get or set the Define components of vector d for AOPT = 2."""  # nopep8
        return self._cards[4].get_value("d3")

    @d3.setter
    def d3(self, value: float) -> None:
        self._cards[4].set_value("d3", value)
