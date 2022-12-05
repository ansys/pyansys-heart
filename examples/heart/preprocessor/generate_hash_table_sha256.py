from ansys.heart.misc.downloader import download_case
import os
import json
import hashlib

if __name__ == "__main__":
    """Creates a reference hash table used to check integrity of data files.

    Note
    ----
    Use this script to download all .tar.gz files from both the Strocchi2020 and
    Rodero2021 public datasets and generate a hash table which is used to
    verify file integrity.
    """
    databases = {"Strocchi2020": 24, "Rodero2021": 20}
    base_folder = "D:\\development\\pyheart-lib\\pyheart-lib\\downloads"
    path_to_hash_table = os.path.join(base_folder, "remote_repo_hash_table_sha256.json")
    if os.path.isfile(path_to_hash_table):
        # read
        fid = open(path_to_hash_table, "r")
        sha256_table: dict = json.load(fid)
        fid.close()
        # convert strings to ints
        for key in sha256_table.keys():
            sha256_table[key] = {int(k): v for k, v in sha256_table[key].items()}
    else:
        sha256_table = {}

    for database in databases:
        num_cases = databases[database]
        download_folder = os.path.join(base_folder, database)
        try:
            sha256_table[database]
        except (KeyError):
            sha256_table[database] = {}

        for ii in range(1, num_cases + 1):
            try:
                sha256_table[database][ii]
                key_exists = True
                print("{0}:{1} Entry already exists".format(database, ii))
                continue
            except (KeyError):
                key_exists = False

            print("\nDownloading... %d" % ii)
            path_to_case = download_case(database, ii, download_folder)
            sha256 = hashlib.sha256(open(path_to_case, "rb").read()).hexdigest()
            sha256_table[database][ii] = sha256

            with open(path_to_hash_table, "w") as outfile:
                json.dump(sha256_table, outfile, indent=4, ensure_ascii=True)
