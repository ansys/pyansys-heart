import typing
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
