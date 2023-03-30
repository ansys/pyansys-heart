import typing
from ansys.dyna.keywords.lib.card import Card, Field
from ansys.dyna.keywords.lib.keyword_base import KeywordBase

class EmEpIsoch(KeywordBase):
    """DYNA EM_EP_ISOCH keyword"""

    keyword = "EM"
    subkeyword = "EP_ISOCH"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field(
                        "idisoch",
                        int,
                        0,
                        10,
                        kwargs.get("idisoch")
                    ),
                    Field(
                        "idepol",
                        int,
                        10,
                        70,
                        kwargs.get("idepol", 0)
                    ),
                    Field(
                        "dplthr",
                        float,
                        20,
                        70,
                        kwargs.get("dplthr", 0)
                    ),
                    Field(
                        "irepol",
                        int,
                        30,
                        70,
                        kwargs.get("irepol", 0)
                    ),
                    Field(
                        "rplthr",
                        float,
                        40,
                        70,
                        kwargs.get("rplthr", 0)
                    ),
                ],
            ),
        ]

    @property
    def idisoch(self) -> typing.Optional[int]:
        """Get or set the ID of the isochrone.
        """ # nopep8
        return self._cards[0].get_value("idisoch")

    @idisoch.setter
    def idisoch(self, value: int) -> None:
        self._cards[0].set_value("idisoch", value)

    @property
    def idepol(self) -> int:
        """Get or set the Flag to activate the computation of depolarization:
        EQ.0: OFF
        EQ.1:ON
        """ # nopep8
        return self._cards[0].get_value("idepol")

    @idepol.setter
    def idepol(self, value: int) -> None:
        self._cards[0].set_value("idepol", value)

    @property
    def dplthr(self) -> float:
        """Get or set the Amplitude threshold used for measuring depolarization.
        """ # nopep8
        return self._cards[0].get_value("dplthr")

    @dplthr.setter
    def dplthr(self, value: float) -> None:
        self._cards[0].set_value("dplthr", value)

    @property
    def irepol(self) -> int:
        """Get or set the Flag to activate the computation of repolarization:
        EQ.0: OFF
        EQ.1:ON
        """ # nopep8
        return self._cards[0].get_value("irepol")

    @irepol.setter
    def irepol(self, value: int) -> None:
        self._cards[0].set_value("irepol", value)

    @property
    def rplthr(self) -> float:
        """Get or set the Amplitude threshold used for measuring repolarization.
        """ # nopep8
        return self._cards[0].get_value("rplthr")

    @rplthr.setter
    def rplthr(self, value: float) -> None:
        self._cards[0].set_value("rplthr", value)

