import pytest

import pandas as pd
import numpy as np

from io import StringIO

from epic.scripts.effective_genome_size import effective_genome_size


@pytest.fixture
def fasta(tmpdir):
    f = tmpdir.join("fasta.file")
    contents = """>chr1
TCAGCAACTGGAAAGGAAACTTTATGTACTGAGTGCTCAGAGTTGTATTA
ACTTTTTTTTTTTTTtgagcagcagcaagatttattgtgaagagtgaaag
aacaaagcttccacagtgtggaaggggacccgagcggtttgccCAGTTGT
ATTAACTTCTAATTCAACACTTTAAGATTCTTAGCATTATTGCAGACAAC
ATcagcttcacaagtgtgtgtcctgtgcagttgaacaagatcccacactt
CCCAGGAGAAACAAGCTTGGCCACTATGCTATCATCAAGTTTCCGCTGAC
>chr2
CACTGAGTCGGCCGGAAGAAGATAGAAGAAAACAACACGCTTGTGTTCAC
ATcagcttcacaagtgtgtgtcctgtgcagttgaacaagatcccacactt
TGTGGATGTTAAAGCCAACAAGCACCAGATCAGACAGGCTGTGAAGAAGC
TCTATGACAGTGATGTGGCCAAGGTCACCACCCTGATTTGTCCTGATAAA
GAGAACAAGGCATATGTTCGACTTGCTCCTGATTATGATGCTTTCGATGT
TGTAACAAAATTGGGATCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"""

    f.write(contents)

    return str(f)


@pytest.fixture
def expected_result():
    pass


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


@pytest.mark.scripts
def test_effective_genome_size(fasta):

    result, result_with_n = effective_genome_size(fasta)
    print(result)
    print(result_with_n)
    assert 0
    # assert isclose(0.9548133595284872, result)
