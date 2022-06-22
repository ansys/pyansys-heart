"""Auto downloads cases from the remote repositories of Strocchi et al 2020, and Cristobal et al 2021
"""
from ansys.heart.custom_logging import logger
import os
from pathlib import Path, PurePath
import requests
import wget
from bs4 import BeautifulSoup
from tqdm import tqdm

URLS = {
    "Strocchi2020" :
    {
        "url": "https://zenodo.org/record/3890034", "num_cases" : 24
    },
    "Cristobal2021" :
    {
        "url" : "https://zenodo.org/record/4590294", "num_cases" : 20
    }
}
DOWNLOAD_DIR = PurePath.joinpath ( Path( __file__ ).parents[3], "downloads" )

def download_cases():
    overwrite_previous = False
    for database_name, subdict in URLS.items():        
        url = subdict["url"]
        num_cases = subdict["num_cases"]

        download_dir = PurePath.joinpath( DOWNLOAD_DIR, database_name )
        if not os.path.isdir(download_dir):
            os.makedirs( download_dir )
        for ii in range(1, num_cases):
            print( "\nDownloading case: {:02d} ...".format(ii) )
            save_path = os.path.join( download_dir, "{:02d}.tar.gz".format(ii) )
            if os.path.isfile(save_path):
                print("File already exists, skipping")
                continue
            download_url = "{:}/files/{:02d}.tar.gz?download=1".format( url, ii )
            wget.download(download_url, save_path)        
    
    return

def unzip_cases():
    """Un-tars the downloaded cases
    """
    import tarfile
    import glob as glob
    for database_name, subdict in URLS.items():        
        download_dir = PurePath.joinpath( DOWNLOAD_DIR, database_name )
        files = glob.glob( os.path.join( download_dir, "*.tar.gz") )
        for file in tqdm(files):
            tar_ball = tarfile.open( file )
            tar_ball.extractall(path = download_dir)    
    return

if __name__ == "__main__":
    download_cases()
    unzip_cases()

    print("Protected")