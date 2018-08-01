import numpy as np

from sys import stderr

cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cpdef _find_islands_cython(const long [::1] bins, const long [::1] chip_counts, const long [::1] input_counts, const double [::1] scores, int distance_allowed, double score_threshold, int window_size):

    cdef int i = 0
    cdef int new_bin = 0
    cdef int old_bin = bins[i]
    cdef int dist = 0
    cdef int start = 0
    cdef int end = 0
    cdef int chip = 0
    cdef int background = 0
    cdef int islands_found = 0
    cdef double score = 0

    cdef int length = len(bins)

    starts_arr = np.zeros(length, dtype=np.long)
    ends_arr = np.zeros(length, dtype=np.long)
    chip_counts_arr = np.zeros(length, dtype=np.long)
    input_counts_arr = np.zeros(length, dtype=np.long)
    scores_arr = np.zeros(length, dtype=np.double)

    cdef long[::1] starts_out
    cdef long[::1] ends_out
    cdef long[::1] chip_out
    cdef long[::1] input_out
    cdef double[::1] scores_out

    starts_out = starts_arr
    ends_out = ends_arr
    chip_out = chip_counts_arr
    input_out = input_counts_arr
    scores_out = scores_arr

    for i in range(length):


        new_bin = bins[i]

        # print("---")
        # print("i", i)
        # print("new_bin", new_bin)
        # print("old_bin", old_bin)
        # print("chip", chip)
        # print("background", background)
        # print("score", score)

        dist = new_bin - old_bin
        # print("distance_allowed", distance_allowed)
        # print("dist", dist)

        if dist > distance_allowed:
            # print("checking if island eligible")
            # print("score_threshold", score_threshold)
            if score >= score_threshold:
                starts_out[islands_found] = start
                ends_out[islands_found] = end + window_size - 1
                chip_out[islands_found] = chip
                input_out[islands_found] = background
                scores_out[islands_found] = score
                islands_found += 1

            start = new_bin
            end = new_bin
            score = scores[i]
            chip = chip_counts[i]
            background = input_counts[i]
        else:
            end = new_bin
            score += scores[i]
            chip += chip_counts[i]
            background += input_counts[i]

        old_bin = new_bin

    # in case we have a new island which we have not included
    if start != new_bin:
        if score > score_threshold:
            starts_out[islands_found] = start
            ends_out[islands_found] = end + window_size - 1
            chip_out[islands_found] = chip
            input_out[islands_found] = background
            scores_out[islands_found] = score
            islands_found += 1

    return starts_arr[:islands_found], ends_arr[:islands_found], chip_counts_arr[:islands_found], input_counts_arr[:islands_found], scores_arr[:islands_found]
