import typing
from ansys.dyna.keywords.lib.card import Card, Field
from ansys.dyna.keywords.lib.keyword_base import KeywordBase

class EmEpTentusscherStimulus(KeywordBase):
    """DYNA EM_EP_TENTUSSCHER_STIMULUS keyword"""

    keyword = "EM"
    subkeyword = "EP_TENTUSSCHER_STIMULUS"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field(
                        "stimid",
                        int,
                        0,
                        10,
                        kwargs.get("stimid")
                    ),
                    Field(
                        "settype",
                        int,
                        10,
                        10,
                        kwargs.get("settype")
                    ),
                    Field(
                        "setid",
                        int,
                        20,
                        10,
                        kwargs.get("setid")
                    ),
                ],
            ),
            Card(
                [
                    Field(
                        "stimstrt",
                        float,
                        0,
                        10,
                        kwargs.get("stimstrt")
                    ),
                    Field(
                        "stimt",
                        float,
                        10,
                        10,
                        kwargs.get("stimt", 1000)
                    ),
                    Field(
                        "stimdur",
                        float,
                        20,
                        10,
                        kwargs.get("stimdur", 2)
                    ),
                    Field(
                        "stimamp",
                        float,
                        30,
                        10,
                        kwargs.get("stimamp", 50)
                    ),
                ],
            ),
        ]

    @property
    def stimid(self) -> typing.Optional[int]:
        """Get or set the ID of the stimulation
        """ # nopep8
        return self._cards[0].get_value("stimid")

    @stimid.setter
    def stimid(self, value: int) -> None:
        self._cards[0].set_value("stimid", value)

    @property
    def settype(self) -> typing.Optional[int]:
        """Get or set the Set type: EQ.1: Segment set, EQ.2: Node set
        """ # nopep8
        return self._cards[0].get_value("settype")

    @settype.setter
    def settype(self, value: int) -> None:
        self._cards[0].set_value("settype", value)

    @property
    def setid(self) -> typing.Optional[int]:
        """Get or set the Node set or segment set ID to be stimulated
        """ # nopep8
        return self._cards[0].get_value("setid")

    @setid.setter
    def setid(self, value: int) -> None:
        self._cards[0].set_value("setid", value)

    @property
    def stimstrt(self) -> typing.Optional[float]:
        """Get or set the Starting time of the stimulation
        """ # nopep8
        return self._cards[1].get_value("stimstrt")

    @stimstrt.setter
    def stimstrt(self, value: float) -> None:
        self._cards[1].set_value("stimstrt", value)

    @property
    def stimt(self) -> float:
        """Get or set the Stimulation period
        """ # nopep8
        return self._cards[1].get_value("stimt")

    @stimt.setter
    def stimt(self, value: float) -> None:
        self._cards[1].set_value("stimt", value)

    @property
    def stimdur(self) -> float:
        """Get or set the Stimulation duration
        """ # nopep8
        return self._cards[1].get_value("stimdur")

    @stimdur.setter
    def stimdur(self, value: float) -> None:
        self._cards[1].set_value("stimdur", value)

    @property
    def stimamp(self) -> float:
        """Get or set the Stimulation amplitude
        """ # nopep8
        return self._cards[1].get_value("stimamp")

    @stimamp.setter
    def stimamp(self, value: float) -> None:
        self._cards[1].set_value("stimamp", value)

