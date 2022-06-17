import os
import shutil 

def clean_directory(directory: str):
    """Cleans the directory by removing it and re-creating it 
    """

    if os.path.isdir(directory):
        shutil.rmtree(directory)
        os.mkdir(directory)
    else:
        os.mkdir(directory)

    return