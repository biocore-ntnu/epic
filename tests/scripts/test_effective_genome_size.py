import pytest

import pandas as pd
import numpy as np

from pyfaidx import Fasta

from io import StringIO

from epic.scripts.effective_genome_size import (effective_genome_size)


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


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


@pytest.mark.scripts
def test_effective_genome_size(fasta):

    result, result_with_n = effective_genome_size(fasta, read_length=30)
    print(result)
    print(result_with_n)
    # assert isclose(0.9548133595284872, result)


"""Not bothering to implement chromosome chunks yet."""

# @pytest.fixture
# def expected_result_chromosome_chunks():
#     return ["TA", "AC", "CG", "GC", "CC", "CC", "CC"]

# @pytest.fixture
# def fasta_simple(tmpdir):

#     f = tmpdir.join("fasta.simple")
#     contents = """>chr1
# TACG
# CCCC"""

#     f.write(contents)

#     return Fasta(str(f))

# @pytest.fixture
# def fasta_pyfaidx(fasta, expected_result_chromosome_chunks):

#     return Fasta(fasta, read_ahead=int(10e3))

# @pytest.mark.scripts
# def test_chromosome_chunks(fasta_simple, expected_result_chromosome_chunks):

#     chromosome = fasta_simple["chr1"]

#     result = list(chromosome_chunks(chromosome, 4, 2))

#     print(result)

#     # assert 0
#     assert result == expected_result_chromosome_chunks

# @pytest.fixture
# def fasta_simple2(tmpdir):

#     f = tmpdir.join("fasta.simple2")
#     contents = """>chr1
# TACG
# CCCC
# GGG"""

#     f.write(contents)

#     return Fasta(str(f))

# @pytest.fixture
# def expected_result_chromosome_chunks2():
#     return ["TA", "AC", "CG", "GC", "CC", "CC", "CC", "CG", "GG", "GG"]

# @pytest.mark.scripts
# def test_chromosome_chunks2(fasta_simple2, expected_result_chromosome_chunks2):

#     chromosome = fasta_simple2["chr1"]

#     result = list(chromosome_chunks(chromosome, 4, 2))

#     print(result)

#     assert result == expected_result_chromosome_chunks2

# @pytest.fixture
# def fasta_simple3(tmpdir):

#     f = tmpdir.join("fasta.simple3")
#     contents = """>chr1
# TACG
# CCCC
# GGG"""

#     f.write(contents)

#     return Fasta(str(f))

# # @pytest.fixture
# # def expected_result_chromosome_chunks3():
# #     return ["TA", "AC", "CG", "GC", "CC", "CC", "CC", "CG", "GG", "GG"]

# # @pytest.mark.scripts
# # def test_chromosome_chunks3(fasta_simple3, expected_result_chromosome_chunks3):

# #     chromosome = fasta_simple3["chr1"]

# #     result = list(chromosome_chunks(chromosome, 3, 2))

# #     print(result)

# #     assert result == expected_result_chromosome_chunks3

# # @pytest.fixture
# # def expected_result():
# #     pass
