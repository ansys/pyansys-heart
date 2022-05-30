import typing
from ansys.dyna.keywords.card import Card, Field
from ansys.dyna.keywords.duplicate_card import DuplicateCard
from ansys.dyna.keywords.keyword_base import KeywordBase

"""
This files contains the keywords that is not supported by dynalib
"""


class SetSegmentTitle(KeywordBase):
    """DYNA SET_SEGMENT keyword"""

    keyword = "SET"
    subkeyword = "SEGMENT_TITLE"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [Field("title", str, 0, 80, kwargs.get("title", "my-title"))],
            ),
            Card(
                [
                    Field("sid", int, 0, 10, kwargs.get("sid")),
                    Field("da1", float, 10, 10, kwargs.get("da1", 0.0)),
                    Field("da2", float, 20, 10, kwargs.get("da2", 0.0)),
                    Field("da3", float, 30, 10, kwargs.get("da3", 0.0)),
                    Field("da4", float, 40, 10, kwargs.get("da4", 0.0)),
                    Field("solver", str, 50, 10, kwargs.get("solver", "MECH")),
                    Field("its", int, 60, 10, kwargs.get("its")),
                ],
            ),
            DuplicateCard(
                [
                    Field("n1", int, 0, 10),
                    Field("n2", int, 10, 10),
                    Field("n3", int, 20, 10),
                    Field("n4", int, 30, 10),
                    Field("a1", float, 40, 10),
                    Field("a2", float, 50, 10),
                    Field("a3", float, 60, 10),
                    Field("a4", float, 70, 10),
                ],
                None,
            ),
        ]

    def _get_title(self):
        return f"*SET_SEGMENT_TITLE"

    @property
    def sid(self) -> typing.Optional[int]:
        """Get or set the Segment set ID. All segment sets should have a unique set ID."""  # nopep8
        return self._cards[1].get_value("sid")

    @sid.setter
    def sid(self, value: int) -> None:
        self._cards[1].set_value("sid", value)

    @property
    def da1(self) -> float:
        """Get or set the First segment attribute default value is 0.0."""  # nopep8
        return self._cards[1].get_value("da1")

    @da1.setter
    def da1(self, value: float) -> None:
        self._cards[1].set_value("da1", value)

    @property
    def da2(self) -> float:
        """Get or set the Second segment attribute default value is 0.0."""  # nopep8
        return self._cards[1].get_value("da2")

    @da2.setter
    def da2(self, value: float) -> None:
        self._cards[1].set_value("da2", value)

    @property
    def da3(self) -> float:
        """Get or set the Third segment attribute default value is 0.0."""  # nopep8
        return self._cards[1].get_value("da3")

    @da3.setter
    def da3(self, value: float) -> None:
        self._cards[1].set_value("da3", value)

    @property
    def da4(self) -> float:
        """Get or set the Fourth segment attribute default value is 0.0."""  # nopep8
        return self._cards[1].get_value("da4")

    @da4.setter
    def da4(self, value: float) -> None:
        self._cards[1].set_value("da4", value)

    @property
    def solver(self) -> str:
        """Get or set the EQ.MECH: mechanics.
         EQ.CESE: CE/SE compressible fluid flow solver.
        EQ.ICFD: Incompressible fluid flow solver.
        """  # nopep8
        return self._cards[1].get_value("solver")

    @solver.setter
    def solver(self, value: str) -> None:
        if value not in ["MECH", "CESE", "ICFD"]:
            raise Exception("""solver must be one of {"MECH","CESE","ICFD"}""")
        self._cards[1].set_value("solver", value)

    @property
    def its(self) -> typing.Optional[int]:
        """Get or set the Define coupling type across different scales in two-scale co-simulation. This flag should only be included for segment sets that provide coupling information in the input file referred to by *INCLUDE_COSIM;
        EQ.1:	Tied contact coupling
        EQ.2 : Solid - in - shell immersed coupling
        """  # nopep8
        return self._cards[1].get_value("its")

    @its.setter
    def its(self, value: int) -> None:
        self._cards[1].set_value("its", value)

    @property
    def segments(self):
        """Gets the table of segments"""
        return self._cards[2].table

    @segments.setter
    def segments(self, df):
        """sets segments from the dataframe df"""
        self._cards[2].table = df


class Mat077H(KeywordBase):
    """DYNA MAT_077_H keyword
        Replace the bug in current version of dynalib

    """

    keyword = "MAT"
    subkeyword = "077_H"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cards = [
            Card(
                [
                    Field("mid", int, 0, 10, kwargs.get("mid")),
                    Field("ro", float, 10, 10, kwargs.get("ro")),
                    Field("pr", float, 20, 10, kwargs.get("pr")),
                    Field("n", int, 30, 10, kwargs.get("n", 0)),
                    Field("nv", int, 40, 10, kwargs.get("nv")),
                    Field("g", float, 50, 10, kwargs.get("g")),
                    Field("sigf", float, 60, 10, kwargs.get("sigf")),
                    Field("ref", float, 70, 10, kwargs.get("ref", 0.0)),
                ],
            ),
            Card(
                [
                    Field("c10", float, 0, 10, kwargs.get("c10")),
                    Field("c01", float, 10, 10, kwargs.get("c01")),
                    Field("c11", float, 20, 10, kwargs.get("c11")),
                    Field("c20", float, 30, 10, kwargs.get("c20")),
                    Field("c02", float, 40, 10, kwargs.get("c02")),
                    Field("c30", float, 50, 10, kwargs.get("c30")),
                    Field("therml", float, 60, 10, kwargs.get("therml")),
                ],
            ),
        ]

    def _get_title(self):
        return f"*MAT_077_H"


if __name__ == "__main__":
    # test isotropic material
    material_iso_kw = Mat077H(
        mid=1,
        ro=1e-6,
        pr=0.499,
        n=0,
        c10=7.45,
    )
    print(material_iso_kw)
