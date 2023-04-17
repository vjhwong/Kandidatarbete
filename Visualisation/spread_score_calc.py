import pandas as pd
import json
from data_extraction_functions import (
    extract_chosen_isolates,
    extract_SIR,
    filter_mic_values,
)
from spread_list_functions import (
    score_mic_spread_list,
    fill_mic_spread_list,
)


def create_mic_spread_dict(
    antibiotic_ranges: dict,
    total_concentration_range: list,
    concentration_to_index_convert: dict,
) -> dict:

    # Dictionary that will contain the mic_spread of each antibiotic
    mic_spread_dict = {
        antibiotic: [0 for _ in range(len(total_concentration_range))]
        for antibiotic in antibiotic_ranges
    }
    
    # Borde kunna förbättras. Kanske combinera värdena till en tuple i samma dict.
    # Blir kanske lite otydligare variabelnamn dock.
    for (_, ranges), (_, mic_spread) in zip(
        antibiotic_ranges.items(), mic_spread_dict.items()
    ):
        # Store lower and upper concentration range of the antibiotic
        lower_limit, upper_limit = ranges["Lower"], ranges["Upper"]

        # Dont leave the mic_spread list untouched for abx with on-scale values
        # on entire range
        if lower_limit == "Min_C" and upper_limit == "Max_C":
            pass
        else:
            # Convert limits from concentration to respective index position
            lower_limit_index = concentration_to_index_convert[lower_limit]
            upper_limit_index = concentration_to_index_convert[upper_limit]
            for index in range(len(mic_spread)):
                if index < lower_limit_index:
                    mic_spread[index] = None
                elif index > upper_limit_index:
                    mic_spread[index] = None
    return mic_spread_dict


def fill_mic_spread_dict(
    antibiotics: list,
    filtered_chosen_isolates_SIR: dict,
    mic_spread_dict: dict,
) -> None:

    for antibiotic in antibiotics:
        mic_spread_list = mic_spread_dict[antibiotic]
        # Creates a set that is converted to a list and then sorted
        unique_mic_values = sorted(
            list(
                {
                    mic_value
                    for _, mic_value, *_ in filtered_chosen_isolates_SIR[antibiotic]
                }
            )
        )
        fill_mic_spread_list(mic_spread_list, unique_mic_values)
    return mic_spread_dict


def score_mic_spread_dict(mic_spread_dict: dict):
    for antibiotic, spread_list in mic_spread_dict.items():
        mic_spread_dict[antibiotic] = (spread_list, score_mic_spread_list(spread_list))


def score_whole_panel(mic_spread_dict: dict) -> int:
    whole_panel_score = 0
    print(f"{'Antibiotic':<30}|{'Score':<7}| Valid spread list")
    for abx, (spread_list, score) in mic_spread_dict.items():
        valid_spread_list = [i for i in spread_list if i is not None]
        print(f"{abx:29} | {round(score, 2):<5} | {valid_spread_list}  ")
        # print(f"{abx:<10}: {valid_spread_list} | score: {score:.2f}")
        whole_panel_score += score

    whole_panel_score /= len(mic_spread_dict)
    print(f"\nWhole panel score: {whole_panel_score:.2f}")


def main():
    # Load files
    chosen_isolates_list = pd.read_csv("Visualisation/Chosen_isolates_list.csv")
    CIB = pd.ExcelFile("Visualisation/CIB_TF-data_AllIsolates_20230302.xlsx")
    matrix_EU = pd.read_excel(CIB, "matrix EU")
    antibiotics_ranges = json.load(open("Visualisation/abx_ranges.json"))

    # Rename a long name for plotting purposes
    # matrix_EU.rename(
    #     columns={"Trimethoprim-sulfamethoxazole": "Trimeth-sulf"}, inplace=True
    # )

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
        "1.0",
        "2.0",
        "4.0",
        "8.0",
        "16.0",
        "32.0",
        "64.0",
        "128.0",
        "256.0",
        "512.0",
        "1024.0",
        "Max C",
    ]

    # Dictionary to go between concentration and indices
    concentration_to_index_convert = {
        concentration: index
        for index, concentration in enumerate(total_concentration_range)
    }

    mic_spread_dict = create_mic_spread_dict(
        antibiotics_ranges, total_concentration_range, concentration_to_index_convert
    )

    # Select isolates
    chosen_isolates = extract_chosen_isolates(chosen_isolates_list, matrix_EU)

    # List of antibiotic names
    antibiotics = list(chosen_isolates.columns[3:])

    # Extract all SIRs for an antibiotic.
    chosen_isolates_SIR = extract_SIR(chosen_isolates, antibiotics)

    # Remove the tuples that have None in their SIR data
    filtered_chosen_isolates_SIR = filter_mic_values(chosen_isolates_SIR)

    mic_spread_dict = fill_mic_spread_dict(
        antibiotics,
        filtered_chosen_isolates_SIR,
        mic_spread_dict,
    )

    score_mic_spread_dict(mic_spread_dict)

    score_whole_panel(mic_spread_dict)


if __name__ == "__main__":
    main()
