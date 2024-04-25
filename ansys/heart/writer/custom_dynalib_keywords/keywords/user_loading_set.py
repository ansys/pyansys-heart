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

from ansys.dyna.keywords.lib.card import Field
from ansys.dyna.keywords.lib.duplicate_card import DuplicateCard
from ansys.dyna.keywords.lib.keyword_base import KeywordBase


class UserLoadingSet(KeywordBase):
    """DYNA USER_LOADING_SET keyword"""

    keyword = "USER"
    subkeyword = "LOADING_SET"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            DuplicateCard(
                [
                    Field("sid", int, 0, 10, kwargs.get("sid")),
                    Field("ltype", str, 10, 10, kwargs.get("ltype")),
                    Field("lcid", int, 20, 10, kwargs.get("lcid")),
                    Field("cid", int, 30, 10, kwargs.get("cid")),
                    Field("sf1", float, 40, 10, kwargs.get("sf1")),
                    Field("sf2", float, 50, 10, kwargs.get("sf2")),
                    Field("sf3", float, 60, 10, kwargs.get("sf3")),
                    Field("iduls", int, 70, 10, kwargs.get("iduls")),
                ],
                None,
            ),
        ]

    @property
    def load_sets(self):
        """Gets the table of load sets"""
        return self._cards[0].table

    @load_sets.setter
    def load_sets(self, df):
        """sets load_sets from the dataframe df"""
        self._cards[0].table = df
