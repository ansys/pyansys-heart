from ansys.dyna.keywords.lib.card import Card, Field
from ansys.dyna.keywords.lib.keyword_base import KeywordBase


class EmEpCellmodelTomek(KeywordBase):
    """DYNA EM_EP_CELLMODEL_TOMEK keyword"""

    keyword = "EM"
    subkeyword = "EP_CELLMODEL_TOMEK"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field("mid", int, 0, 10, kwargs.get("mid", 1)),
                    Field("phiendmid", float, 10, 10, kwargs.get("phiendmid", 0.17)),
                    Field("phimidepi", float, 20, 10, kwargs.get("phimidepi", 0.58)),
                ],
            ),
        ]

    @property
    def mid(self) -> int:
        """Get or set the Material ID"""  # nopep8
        return self._cards[0].get_value("mid")

    @mid.setter
    def mid(self, value: int) -> None:
        self._cards[0].set_value("mid", value)

    @property
    def phiendmid(self) -> float:
        """Get or set the Phi endocardium > mid"""  # nopep8
        return self._cards[0].get_value("phiendmid")

    @phiendmid.setter
    def phiendmid(self, value: float) -> None:
        self._cards[0].set_value("phiendmid", value)

    @property
    def phimidepi(self) -> float:
        """Get or set the Phi mid > epicardium"""  # nopep8
        return self._cards[0].get_value("phimidepi")

    @phimidepi.setter
    def phimidepi(self, value: float) -> None:
        self._cards[0].set_value("phimidepi", value)
