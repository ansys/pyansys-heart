import typing
from ansys.dyna.keywords.lib.card import Card, Field
from ansys.dyna.keywords.lib.keyword_base import KeywordBase

class UserLoadingSet(KeywordBase):
    """DYNA USER_LOADING_SET keyword"""

    keyword = "USER"
    subkeyword = "LOADING_SET"

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
                    Field(
                        "ltype",
                        str,
                        10,
                        10,
                        kwargs.get("ltype")
                    ),
                    Field(
                        "lcid",
                        int,
                        20,
                        10,
                        kwargs.get("lcid")
                    ),
                    Field(
                        "cid",
                        int,
                        30,
                        10,
                        kwargs.get("cid")
                    ),
                    Field(
                        "sf1",
                        float,
                        40,
                        10,
                        kwargs.get("sf1")
                    ),
                    Field(
                        "sf2",
                        float,
                        50,
                        10,
                        kwargs.get("sf2")
                    ),
                    Field(
                        "sf3",
                        float,
                        60,
                        10,
                        kwargs.get("sf3")
                    ),
                    Field(
                        "iduls",
                        int,
                        70,
                        10,
                        kwargs.get("iduls")
                    ),
                ],
            ),
        ]

    @property
    def sid(self) -> typing.Optional[int]:
        """Get or set the ID of the set to which user-defined loading will be applied. Set type depends on the type of loading, see LTYPE.
        """ # nopep8
        return self._cards[0].get_value("sid")

    @sid.setter
    def sid(self, value: int) -> None:
        self._cards[0].set_value("sid", value)

    @property
    def ltype(self) -> typing.Optional[str]:
        """Get or set the Loading type
        """ # nopep8
        return self._cards[0].get_value("ltype")

    @ltype.setter
    def ltype(self, value: str) -> None:
        self._cards[0].set_value("ltype", value)

    @property
    def lcid(self) -> typing.Optional[int]:
        """Get or set the Load curve, a function of time. Its current value, crv, is passed to user subroutine LOADSETUD.
        """ # nopep8
        return self._cards[0].get_value("lcid")

    @lcid.setter
    def lcid(self, value: int) -> None:
        self._cards[0].set_value("lcid", value)

    @property
    def cid(self) -> typing.Optional[int]:
        """Get or set the Optional coordinate system along which scale factors SFi is defined. Global system is the default system.
        """ # nopep8
        return self._cards[0].get_value("cid")

    @cid.setter
    def cid(self, value: int) -> None:
        self._cards[0].set_value("cid", value)

    @property
    def sf1(self) -> typing.Optional[float]:
        """Get or set the Scale factor of loading magnitude, when LTYPE
        """ # nopep8
        return self._cards[0].get_value("sf1")

    @sf1.setter
    def sf1(self, value: float) -> None:
        self._cards[0].set_value("sf1", value)

    @property
    def sf2(self) -> typing.Optional[float]:
        """Get or set the Scale factor of loading magnitude, when LTYPE
        """ # nopep8
        return self._cards[0].get_value("sf2")

    @sf2.setter
    def sf2(self, value: float) -> None:
        self._cards[0].set_value("sf2", value)

    @property
    def sf3(self) -> typing.Optional[float]:
        """Get or set the Scale factor of loading magnitude, when LTYPE
        """ # nopep8
        return self._cards[0].get_value("sf3")

    @sf3.setter
    def sf3(self, value: float) -> None:
        self._cards[0].set_value("sf3", value)

    @property
    def iduls(self) -> typing.Optional[int]:
        """Get or set the Each USER_LOADING_SET can be assigned a unique ID, which is passed to user subroutine LOADSETUD and allows multiple loading definitions by using a single user subroutine, LOADSETUD. If no value is input, LS-DYNA will assign a sequence number to each USER_LOADING_SET based on its definition sequence.
        """ # nopep8
        return self._cards[0].get_value("iduls")

    @iduls.setter
    def iduls(self, value: int) -> None:
        self._cards[0].set_value("iduls", value)

