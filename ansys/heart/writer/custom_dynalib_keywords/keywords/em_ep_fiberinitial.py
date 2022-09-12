# flake8: noqa: E501
import typing

from ansys.dyna.keywords.lib.card import Card, Field
from ansys.dyna.keywords.lib.keyword_base import KeywordBase


class EmEpFiberinitial(KeywordBase):
    """DYNA EM_EP_FIBERINITIAL keyword"""

    keyword = "EM"
    subkeyword = "EP_FIBERINITIAL"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field("id", int, 0, 10, kwargs.get("id")),
                    Field("partid", int, 10, 10, kwargs.get("partid")),
                    Field("stype", int, 20, 10, kwargs.get("stype", 1)),
                    Field("ssid1", int, 30, 10, kwargs.get("ssid1")),
                    Field("ssid2", int, 40, 10, kwargs.get("ssid2")),
                ],
            ),
        ]

    @property
    def id(self) -> typing.Optional[int]:
        """Get or set the ID of the Laplace system to solve (define new id with each new line)"""  # nopep8
        return self._cards[0].get_value("id")

    @id.setter
    def id(self, value: int) -> None:
        self._cards[0].set_value("id", value)

    @property
    def partid(self) -> typing.Optional[int]:
        """Get or set the Part id on which the system is solved"""  # nopep8
        return self._cards[0].get_value("partid")

    @partid.setter
    def partid(self, value: int) -> None:
        self._cards[0].set_value("partid", value)

    @property
    def stype(self) -> int:
        """Get or set the Segment type: eq = 1: node-set, eq=2: segment-set"""  # nopep8
        return self._cards[0].get_value("stype")

    @stype.setter
    def stype(self, value: int) -> None:
        self._cards[0].set_value("stype", value)

    @property
    def ssid1(self) -> typing.Optional[int]:
        """Get or set the Set on which a potential of value 1 is prescribed"""  # nopep8
        return self._cards[0].get_value("ssid1")

    @ssid1.setter
    def ssid1(self, value: int) -> None:
        self._cards[0].set_value("ssid1", value)

    @property
    def ssid2(self) -> typing.Optional[int]:
        """Get or set the Set on which a potential of value 0 is prescribed"""  # nopep8
        return self._cards[0].get_value("ssid2")

    @ssid2.setter
    def ssid2(self, value: int) -> None:
        self._cards[0].set_value("ssid2", value)
