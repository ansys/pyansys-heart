import typing
from ansys.dyna.keywords.lib.card import Card, Field
from ansys.dyna.keywords.lib.keyword_base import KeywordBase


class EmEpEcg(KeywordBase):
    """DYNA EM_EP_ECG keyword"""

    keyword = "EM"
    subkeyword = "EP_ECG"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field("ecgid", int, 0, 10, kwargs.get("ecgid")),
                    Field("psid", int, 10, 10, kwargs.get("psid")),
                ],
            ),
        ]

    @property
    def ecgid(self) -> typing.Optional[int]:
        """Get or set the ID of the ECG computation
        """  # nopep8
        return self._cards[0].get_value("ecgid")

    @ecgid.setter
    def ecgid(self, value: int) -> None:
        self._cards[0].set_value("ecgid", value)

    @property
    def psid(self) -> typing.Optional[int]:
        """Get or set the Point set ID containing the list of virtual points on which the pseudo-ECGs are computed
        """  # nopep8
        return self._cards[0].get_value("psid")

    @psid.setter
    def psid(self, value: int) -> None:
        self._cards[0].set_value("psid", value)
