import pytest

from io import StringIO

from epic.utils.find_readlength import (find_readlength,
                                        get_closest_readlength)

__author__ = "Endre Bakken Stovner https://github.com/endrebak/"
__license__ = "MIT"


def test_get_closest_readlength():

    assert 36 == get_closest_readlength(37)
    assert 36 == get_closest_readlength(35)
    assert 36 == get_closest_readlength(-1)
    assert 50 == get_closest_readlength(50)
    assert 50 == get_closest_readlength(62)
    assert 100 == get_closest_readlength(1500)


def test_find_readlength(args_200):

    result = find_readlength(args_200)
    assert result == 25
