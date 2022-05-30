import typing
from ansys.dyna.keywords.card import Card, Field
from ansys.dyna.keywords.duplicate_card import DuplicateCard
from ansys.dyna.keywords.keyword_base import KeywordBase

"""
This files contains the keywords that is not supported by dynalib
"""

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
