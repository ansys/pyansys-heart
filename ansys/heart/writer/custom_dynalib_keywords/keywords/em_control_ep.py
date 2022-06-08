import typing
from ansys.dyna.keywords.lib.card import Card, Field
from ansys.dyna.keywords.lib.keyword_base import KeywordBase

class EmControlEp(KeywordBase):
    """DYNA EM_CONTROL_EP keyword"""

    keyword = "EM"
    subkeyword = "CONTROL_EP"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field(
                        "solvetype",
                        int,
                        0,
                        10,
                        kwargs.get("solvetype", 4)
                    ),
                    Field(
                        "numsplit",
                        int,
                        10,
                        10,
                        kwargs.get("numsplit", 5)
                    ),
                    Field(
                        "actusig",
                        int,
                        20,
                        10,
                        kwargs.get("actusig", 100000000)
                    ),
                ],
            ),
        ]

    @property
    def solvetype(self) -> int:
        """Get or set the ?
        """ # nopep8
        return self._cards[0].get_value("solvetype")

    @solvetype.setter
    def solvetype(self, value: int) -> None:
        self._cards[0].set_value("solvetype", value)

    @property
    def numsplit(self) -> int:
        """Get or set the ?
        """ # nopep8
        return self._cards[0].get_value("numsplit")

    @numsplit.setter
    def numsplit(self, value: int) -> None:
        self._cards[0].set_value("numsplit", value)

    @property
    def actusig(self) -> int:
        """Get or set the ?
        """ # nopep8
        return self._cards[0].get_value("actusig")

    @actusig.setter
    def actusig(self, value: int) -> None:
        self._cards[0].set_value("actusig", value)

