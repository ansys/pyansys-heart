import typing
from ansys.dyna.keywords.card import Card, Field
from ansys.dyna.keywords.option import Options, Option
from ansys.dyna.keywords.duplicate_card import DuplicateCard
from ansys.dyna.keywords.keyword_base import KeywordBase

class SetSegmentAdd_custom(KeywordBase):
    """DYNA SET_SEGMENT_ADD keyword"""

    keyword = "SET"
    subkeyword = "SEGMENT_ADD"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field(
                        "sid",
                        int,
                        0,
                        10,
                        kwargs.get("sid")
                    ),
                ],
            ),
            DuplicateCard(
                [
                    Field("ssid1", int, 0, 10),
                    Field("ssid2", int, 10, 10),
                    Field("ssid3", int, 20, 10),
                    Field("ssid4", int, 30, 10),
                    Field("ssid5", int, 40, 10),
                    Field("ssid6", int, 50, 10),
                    Field("ssid7", int, 60, 10),
                    Field("ssid8", int, 70, 10),
                ],
                None,
            ),
        ]
        self._options = Options([
            Option(
                name = "TITLE",
                card_order = -1,
                title_order = 1,
                cards = [
                    Card(
                        [
                            Field(
                                "title",
                                str,
                                0,
                                80,
                                kwargs.get("title")
                            ),
                        ],
                    ),
                ]
            ),
        ])

    @property
    def sid(self) -> typing.Optional[int]:
        """Get or set the Segment set ID. All segment sets should have a unique set ID.
        """ # nopep8
        return self._cards[0].get_value("sid")

    @sid.setter
    def sid(self, value: int) -> None:
        self._cards[0].set_value("sid", value)

    @property
    def segment_ids(self):
        '''Gets the table of segment ids'''
        return self._cards[1].table

    @segment_ids.setter
    def segment_ids(self, df):
        '''sets list of segment ids from the dataframe df'''
        self._cards[1].table = df

    @property
    def title(self) -> typing.Optional[str]:
        """Get or set the Additional title line
        """ # nopep8
        return self._options["TITLE"].cards[0].get_value("title")

    @title.setter
    def title(self, value: str) -> None:
        self._options["TITLE"].cards[0].set_value("title", value)

