import numpy as np


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

    if valid_list[0] == 0:
        edge_penalty += 0.5
    if valid_list[-1] == 0:
        edge_penalty += 0.5

    return edge_penalty


def score_mic_spread_list(fill_list: list) -> float:
    valid_list = [i for i in fill_list if i is not None]

    if len(valid_list) <= 1:
        raise ValueError("Length of valid list must be greater than 1")

    penalty = 0
    penalty += count_gap_length(valid_list)
    penalty += check_edges(valid_list)

    # Lägg till att dela med längden av listan
    return penalty


def fill_mic_spread_list(
    mic_spread_list: list, unique_mic_values: set, concentration_to_index_convert: dict
) -> None:

    for mic_value in unique_mic_values:
        index = int(np.log2(mic_value) + 10)
        mic_spread_list[index] = 1
