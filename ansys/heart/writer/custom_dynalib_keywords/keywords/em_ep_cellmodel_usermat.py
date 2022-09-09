from ansys.dyna.keywords.lib.card import Card, Field
from ansys.dyna.keywords.lib.keyword_base import KeywordBase


class EmEpCellmodelUsermat(KeywordBase):
    """DYNA EM_EP_CELLMODEL_USERMAT keyword"""

    keyword = "EM"
    subkeyword = "EP_CELLMODEL_USERMAT"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field("mid", int, 0, 10, kwargs.get("mid", 1)),
                ],
            ),
        ]

    @property
    def mid(self) -> int:
        """Get or set the Material ID. A unique number must be specified (see *PART)."""  # nopep8
        return self._cards[0].get_value("mid")

    @mid.setter
    def mid(self, value: int) -> None:
        self._cards[0].set_value("mid", value)
