""" Global settings for the logger module. """

import logging
import sys

logging.basicConfig(
    level=logging.DEBUG,
    format=
    '%(message)s (File: %(module)s, Log level: %(levelname)s, Time: %(asctime)s )',
    datefmt='%a, %d %b %Y %H:%M:%S',
    stream=sys.stderr)
