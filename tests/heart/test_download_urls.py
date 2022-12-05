from ansys.heart.misc.downloader import download_case, unzip_case, format_download_urls

from .conftest import get_assets_folder, get_workdir

import pytest
import validators


@pytest.mark.parametrize(
    "database_name",
    ["Strocchi2020", "Rodero2021"],
)
def test_download_urls(database_name):
    """Test if URL still valid and exists."""
    all_download_urls = format_download_urls()
    for case_num in all_download_urls[database_name]:
        url = all_download_urls[database_name][case_num]
        assert validators.url(url), "No valid URL for Case {0} of database {1}".format(
            case_num, database_name
        )


# how to check whether urls still contain the right data?
