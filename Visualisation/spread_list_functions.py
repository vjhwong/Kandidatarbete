import numpy as np

"""
Functions used to calculate spread of bacterial isolates MIC values for an antibiotic using method developed by
Victor Wong and Hanna Malmvall 12 April 20223. Also contains one function that fills a
spread list given the unique MIC values of an antibiotic.  

Rules:
1. Count number of empty values in a row and subtract one.
2. Add 0.5 penalty for each empty edge. 

Author: Victor Wong
Date: 16 April 2023
"""


def count_gap_length(valid_list: list) -> float:
    total_gap_length = 0
    gap_length = 0

    for value in valid_list:
        if value == 0:
            gap_length += 1
        else:  # if current value == 1
            if gap_length >= 2:
                total_gap_length += gap_length - 1
            gap_length = 0
    # If the last values were also 0s and the gap length is at least 2.
    # Add the length of the gap to the total gap length
    if gap_length >= 2:
        total_gap_length += gap_length - 1
    return total_gap_length


def check_edges(valid_list: list) -> float:
    edge_penalty = 0

    # Add 0.5 penalty if the edge values are empty.
    if valid_list[0] == 0:  # Check first value
        edge_penalty += 0.5
    if valid_list[-1] == 0:  # Check last value
        edge_penalty += 0.5

    return edge_penalty


def score_mic_spread_list(fill_list: list) -> float:
    valid_list = [i for i in fill_list if i is not None]

    if len(valid_list) <= 1:
        raise ValueError("Length of valid list must be greater than 1")

    penalty = 0
    penalty += count_gap_length(valid_list)
    penalty += check_edges(valid_list)
    print(penalty)
    score = 1 - penalty / len(valid_list)
    # Lägg till att dela med längden av listan
    return score


def fill_mic_spread_list(
    mic_spread_list: list,
    unique_mic_values: list,
) -> None:

    for mic_value in unique_mic_values:
        index = int(np.log2(mic_value) + 10)
        mic_spread_list[index] = 1
