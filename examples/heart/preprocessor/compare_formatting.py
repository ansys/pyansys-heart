"""Example to pre-process data from Strocchi2020 and Cristobal2021"""
import os
import pathlib

if __name__ == "__main__":

    """Full Heart example

    1. Extracts simulation mesh
    2. Writes files for mechanics, zero-pressure, fiber generation, and purkinje
    """

    path_to_case = "D:\\development\\pyheart-lib\\pyheart-lib\\downloads\\Strocchi2020\\02\\02.case"
    workdir_bad = os.path.join(
        pathlib.Path(path_to_case).parent, "BiVentricleRefactoredBadFormatting_2mm"
    )

    workdir_good = os.path.join(
        pathlib.Path(path_to_case).parent, "BiVentricleRefactoredGoodFormatting_2mm"
    )
    path_to_model = os.path.join(workdir_good, "heart_model.pickle")

    from ansys.dyna.keywords import Deck, keywords
    import numpy as np

    keywords.SectionShell()

    deck_good = Deck()
    deck_bad = Deck()
    files_to_compare = ["mechanics\\pericardium.k", "mechanics\\boundary_conditions.k"]
    for file in files_to_compare:
        deck_good.import_file(os.path.join(workdir_good, file))
        deck_bad.import_file(os.path.join(workdir_bad, file))
        if "pericardium" in file:
            deck_bad.keywords[4].elements["s"]
            deck_good.keywords[4].elements["s"]
        if "boundary_condition" in file:
            for ii, kw in enumerate(deck_good.keywords):
                print(kw.keyword + "_" + kw.subkeyword + " {}".format(ii))
                if "SD_ORIENTATION" in kw.subkeyword:

                    print(
                        np.max(
                            np.abs(deck_good.keywords[ii].vectors - deck_bad.keywords[ii].vectors)
                        )
                    )
                elif "ELEMENT" in kw.keyword and "DISCRETE" in kw.subkeyword:
                    print(
                        np.max(
                            np.abs(deck_good.keywords[ii].elements - deck_bad.keywords[ii].elements)
                        )
                    )

    pass
