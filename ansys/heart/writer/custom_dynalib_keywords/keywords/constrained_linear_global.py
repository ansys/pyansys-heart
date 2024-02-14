import typing
from ansys.dyna.keywords.lib.card import Card, Field
from ansys.dyna.keywords.lib.duplicate_card import DuplicateCard
from ansys.dyna.keywords.lib.keyword_base import KeywordBase
from ansys.dyna.keywords.lib.option import Option, Options

class ConstrainedLinearGlobal(KeywordBase):
    """DYNA CONSTRAINED_LINEAR_GLOBAL keyword"""

    keyword = "CONSTRAINED"
    subkeyword = "LINEAR_GLOBAL"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field(
                        "licd",
                        int,
                        0,
                        10,
                        kwargs.get("licd")
                    ),
                ],
            ),
            DuplicateCard(
                [
                    Field("nid", int, 0, 10, kwargs.get("nid")),
                    Field("dof", int, 10, 10, kwargs.get("dof")),
                    Field("coef", float, 20, 10, kwargs.get("coef")),
                ],
                None,
            ),
        ]

    @property
    def linear_constraints(self):
        """Gets the table of linear constraints"""
        return self._cards[1].table

    @linear_constraints.setter
    def linear_constraints(self, df):
        """sets linear constraints from the dataframe df"""
        self._cards[1].table = df


    @property
    def licd(self) -> typing.Optional[int]:
        """Get or set the Linear constraint definition ID. This ID can be used to identify a set to which this constraint is a member.
        """ # nopep8
        return self._cards[0].get_value("licd")

    @licd.setter
    def licd(self, value: int) -> None:
        self._cards[0].set_value("licd", value)
