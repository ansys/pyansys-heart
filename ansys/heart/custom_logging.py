import logging
from time import strftime

strftime("%H")

# create logger
logger = logging.getLogger("simple_example")
logger.setLevel(logging.DEBUG)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter
# dateformat = '%b %d %Y %H:%M:%S'
dateformat = "%Y/%m/%d %H:%M:%S"
formatter = logging.Formatter(
    "%(asctime)s - %(levelname)s - %(message)s", datefmt=dateformat
)

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)


if __name__ == "__main__":
    print()
    # log = CustomLogger( to_stdout = True )
