import pandas as pd
import numpy as np
import json
from pprint import pprint
from data_extraction_functions import (
    extract_chosen_isolates,
    extract_mic_data,
    extract_SIR,
    filter_mic_values,
)
from gap_count_functions import count_gap_length, check_edges, score_fill_list


def create_spread_dict(
    antibiotic_ranges: dict, total_concentration_range: list
) -> dict:

    # Dictionary that will contain the spread of each antibiotic
    spread_dict = {
        antibiotic: [0 for _ in range(len(total_concentration_range))]
        for antibiotic in antibiotic_ranges
    }

    # Dictionary to go between concentration and indices
    concentration_to_index_convert = {
        concentration: index
        for index, concentration in enumerate(total_concentration_range)
    }

    for (_, ranges), (_, spread) in zip(antibiotic_ranges.items(), spread_dict.items()):
        # Store lower and upper concentration range of the antibiotic
        lower_limit, upper_limit = ranges["Lower"], ranges["Upper"]

        # Dont leave the spread list untouched for abx with on-scale values
        # on entire range
        if lower_limit == "Min_C" and upper_limit == "Max_C":
            pass
        else:
            # Convert limits from concentration to respective index position
            lower_limit_index = concentration_to_index_convert[lower_limit]
            upper_limit_index = concentration_to_index_convert[upper_limit]
            for index in range(len(spread)):
                if index < lower_limit_index:
                    spread[index] = None
                elif index > upper_limit_index:
                    spread[index] = None
    return spread_dict


def main():
    # Load files
    chosen_isolates_list = pd.read_csv("Visualisation/Chosen_isolates_list.csv")
    CIB = pd.ExcelFile("Visualisation/CIB_TF-data_AllIsolates_20230302.xlsx")
    matrix_EU = pd.read_excel(CIB, "matrix EU")
    antibiotics_ranges = json.load(open("Visualisation/abx_ranges.json"))

    # Rename a long name for plotting purposes
    matrix_EU.rename(
        columns={"Trimethoprim-sulfamethoxazole": "Trimeth-sulf"}, inplace=True
    )

    total_concentration_range = [
        "Min C",
        "0.00195",
        "0.00391",
        "0.00781",
        "0.01563",
        "0.03125",
        "0.0625",
        "0.125",
        "0.25",
        "0.5",
        "1",
        "2",
        "4",
        "8",
        "16",
        "32",
        "64",
        "128",
        "256",
        "512",
        "1024",
        "Max C",
    ]

    spread_dict = create_spread_dict(antibiotics_ranges, total_concentration_range)

    # # Select isolates
    # chosen_isolates = extract_chosen_isolates(chosen_isolates_list, matrix_EU)

    # # List of antiiotic names
    # antibiotics = list(chosen_isolates.columns[3:])

    # # Extract all SIRs for an antibiotic.
    # chosen_isolates_SIR = extract_SIR(chosen_isolates, antibiotics)

    # # Remove the tuples that have None in their SIR data
    # filtered_chosen_isolates_SIR = filter_mic_values(chosen_isolates_SIR)

    # # Extract the mic-values of each isolate for each antibiotic.
    # mic_data = extract_mic_data(filtered_chosen_isolates_SIR, antibiotics)


if __name__ == "__main__":
    main()
