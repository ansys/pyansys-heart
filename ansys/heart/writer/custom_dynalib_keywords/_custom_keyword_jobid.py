import typing
from ansys.dyna.keywords.card import Card, Field
from ansys.dyna.keywords.keyword_base import KeywordBase

class KeywordJobid(KeywordBase):
    """DYNA KEYWORD_KEYWORD_JOBID keyword"""

    keyword = "KEYWORD"
    subkeyword = "JOBID"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field(
                        "jobname",
                        str,
                        0,
                        256,
                        kwargs.get("jobname")
                    ),
                ],
            ),
        ]

    @property
    def jobname(self) -> typing.Optional[str]:
        """Get or set the job name
        """ # nopep8
        return self._cards[1].get_value("jobname")

    @jobname.setter
    def jobname(self, value: str) -> None:
        self._cards[1].set_value("jobname", value)

