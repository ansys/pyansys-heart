import typing
from ansys.dyna.keywords.lib.card import Card, Field
from ansys.dyna.keywords.lib.keyword_base import KeywordBase

class EmEpTentusscherStimulus2(KeywordBase):
    """DYNA EM_EP_TENTUSSCHER_STIMULUS2 keyword"""

    keyword = "EM"
    subkeyword = "EP_TENTUSSCHER_STIMULUS2"

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
                    Field(
                        "lcid",
                        int,
                        30,
                        10,
                        kwargs.get("lcid")
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
    def lcid(self) -> typing.Optional[int]:
        """Get or set the load curve to use for stimulation, where the first coordinate represents time and the second represents the stim. amplitude
        """ # nopep8
        return self._cards[0].get_value("lcid")

    @lcid.setter
    def lcid(self, value: int) -> None:
        self._cards[0].set_value("lcid", value)

