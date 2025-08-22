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

# from ansys.dyna.core.lib.duplicate_card import DuplicateCard
from ansys.dyna.core.lib.keyword_base import KeywordBase

"""
This files contains the keywords that is not supported by the PyDYNA keywords module
"""


class Mat077H(KeywordBase):
    """LS-DYNA MAT_077_H keyword
    Replace the bug in current version of PyDYNA keywords module

    """

    keyword = "MAT"
    subkeyword = "077_H"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field("mid", int, 0, 10, kwargs.get("mid")),
                    Field("ro", float, 10, 10, kwargs.get("ro")),
                    Field("pr", float, 20, 10, kwargs.get("pr")),
                    Field("n", int, 30, 10, kwargs.get("n", 0)),
                    Field("nv", int, 40, 10, kwargs.get("nv")),
                    Field("g", float, 50, 10, kwargs.get("g")),
                    Field("sigf", float, 60, 10, kwargs.get("sigf")),
                    Field("ref", float, 70, 10, kwargs.get("ref", 0.0)),
                ],
            ),
            Card(
                [
                    Field("c10", float, 0, 10, kwargs.get("c10")),
                    Field("c01", float, 10, 10, kwargs.get("c01")),
                    Field("c11", float, 20, 10, kwargs.get("c11")),
                    Field("c20", float, 30, 10, kwargs.get("c20")),
                    Field("c02", float, 40, 10, kwargs.get("c02")),
                    Field("c30", float, 50, 10, kwargs.get("c30")),
                    Field("therml", float, 60, 10, kwargs.get("therml")),
                ],
            ),
        ]

    def _get_title(self):
        return f"*MAT_077_H"

    @property
    def mid(self) -> typing.Optional[int]:
        """Get or set the Material ID: refers to MID in the *PART card."""  # nopep8
        return self._cards[0].get_value("mid")

    @mid.setter
    def mid(self, value: int) -> None:
        self._cards[0].set_value("mid", value)


if __name__ == "__main__":
    # test isotropic material
    material_iso_kw = Mat077H(
        mid=1,
        ro=1e-6,
        pr=0.499,
        n=0,
        c10=7.45,
    )
    print(material_iso_kw)
