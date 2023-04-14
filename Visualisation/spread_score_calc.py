import pandas as pd
import numpy as np
import json
from pprint import pprint


def extract_chosen_isolates(
    chosen_isolates: pd.DataFrame, matrix_EU: pd.DataFrame
) -> pd.DataFrame:
    """
    Select the chosen isolates from the script. Return a DataFrame only containing the rows of selected isolates
    """
    chosen_rows = matrix_EU["Isolate"].isin(chosen_isolates["Isolate"])
    return matrix_EU[chosen_rows]


def find_digits(SIR: str) -> int:
    """Find numbers in a string"""
    digit = ""
    for character in SIR:
        if character.isdigit() or character == ".":
            digit += character
    return float(digit)


def get_scale(SIR: str) -> bool:
    """Get on or off scale. True == on-scale"""
    if "=" in SIR:
        return True
    elif "<" in SIR or ">" in SIR:
        return False
    else:
        raise ValueError("Not a valid SIR")


def parse_SIR(SIR: str) -> bool:
    """
    Find the isolates with valid SIRs. Not 'Missing BP'
    and not 'nip'.
    """
    if SIR.startswith("Missing BP"):
        return False
    if SIR == "nip":
        return False
    return True


def extract_SIR(chosen_isolates: pd.DataFrame, antibiotics: list) -> dict:
    """
    Extract all SIRs for an antibiotic. Returns a dictionary
    with antibiotcs as keys and lists of the isolates and their
    SIRs in tuples as value.
    """
    chosen_isolates_SIR = {antibiotic: [] for antibiotic in antibiotics}

    for index, row in chosen_isolates.iterrows():
        isolate, pathogen, antibiotic_SIR = row[0], row[1], list(row[3:].items())
        for antibiotic, SIR in antibiotic_SIR:
            if parse_SIR(SIR):
                mic_category = SIR[0]
                mic = find_digits(SIR)
                scale = get_scale(SIR)
                chosen_isolates_SIR[antibiotic].append(
                    (isolate, mic, mic_category, scale, pathogen)
                )
            else:
                # If SIR = "Missing BP" or "nip"
                chosen_isolates_SIR[antibiotic].append(
                    (isolate, SIR, None, None, pathogen)
                )
    return chosen_isolates_SIR


def filter_mic_values(chosen_isolates_SIR: dict) -> None:
    """
    Remove the tuples that have None in their SIR data
    """
    for antibiotic, SIR_data in chosen_isolates_SIR.items():
        # tup = (isolate, mic_value, mic_category, scale, pathogen)
        chosen_isolates_SIR[antibiotic] = list(
            (tup for tup in SIR_data if tup[2] is not None)
        )
    return chosen_isolates_SIR


def extract_mic_data(chosen_isolates_SIR: dict, antibiotics: list) -> list:
    """
    Extract the mic-values of each isolate for each antibiotic.
    Returns a nested list. Each list represents the mic-values of
    all isolates for an antibiotic.
    """
    mic_values = []
    # Iterate over all antibioticsz
    for antibiotic in antibiotics:
        # Create a list to hold the mic-values of isolates for that abx
        antibiotic_mic_values = []
        # Get value of current abx. List of (isolate, mic_value, mic_category)
        SIR_data = chosen_isolates_SIR[antibiotic]
        for isolate, mic_value, mic_category, scale, pathogen in SIR_data:
            antibiotic_mic_values.append(
                (isolate, np.log2(mic_value), mic_category, scale, pathogen)
            )
        mic_values.append(antibiotic_mic_values)
    return mic_values


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
        lower_limit, upper_limit = ranges["Lower"], ranges["Upper"]
        if lower_limit == "Min_C" and upper_limit == "Max_C":
            pass
        else:
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
    print(spread_dict)
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
