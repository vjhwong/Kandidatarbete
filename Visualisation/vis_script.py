import pandas as pd
import numpy as np
import plotly.express as px


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


def create_plot_df(
    antibiotics: list, mic_data: list, x_jitter: float = 0.15, y_jitter: float = 0.05
) -> pd.DataFrame:
    """Create dataframe used for plotting"""

    mic_dict = {
        "S": ("limegreen", "Sensitive"),
        "I": ("gold", "Intermediate"),
        "R": ("tomato", "Resistant"),
    }

    # Set ticks of x axis
    x_axis = [i for i in range(len(antibiotics))]

    # Initialize lists to hold values used for plotting
    x_values, y_values = [], []
    SIR_category_list = []
    isolate_names = []
    on_off_scale = []
    pathogen_list = []

    for x_value, abx_mic_data in zip(x_axis, mic_data):
        for isolate, mic_value, SIR_category, scale, pathogen in abx_mic_data:
            # Add random noise to avoid overlapping
            x_value_jitter = x_value + np.random.uniform(-x_jitter, x_jitter)
            mic_value_jitter = mic_value + np.random.uniform(-y_jitter, y_jitter)
            # Add data to the lists
            isolate_names.append(isolate)
            x_values.append(x_value_jitter)
            SIR_category_list.append(mic_dict[SIR_category][1])
            on_off_scale.append(scale)
            pathogen_list.append(pathogen)
            if scale is True:
                y_values.append(mic_value_jitter)
            elif scale is False:
                if SIR_category == "S":
                    y_values.append(-10)
                elif SIR_category == "R":
                    y_values.append(11)
                else:
                    raise ValueError(
                        f"SIR Category must be either S or R, not: {SIR_category}"
                    )
            else:
                raise ValueError(
                    f"scale must be Boolean value, not {type(scale)}. Current value: {scale}"
                )

    # Create a DF used for plotting
    plot_df = pd.DataFrame(
        {
            "Antibiotics": x_values,
            "Log2(MIC-value)": y_values,
            "Isolate names": isolate_names,
            "SIR": SIR_category_list,
            "Scale": on_off_scale,
            "Pathogen": pathogen_list,
        },
        index=np.arange(len(x_values)),
    )

    return plot_df


def plotly_dotplot(plot_df: pd.DataFrame, antibiotics: list) -> None:

    # Set ticks of x axis
    x_axis = [i for i in range(len(antibiotics))]

    y_axis_ticktext = [
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

    # plot
    fig = px.scatter(
        plot_df,
        x="Antibiotics",
        y="Log2(MIC-value)",
        hover_name="Isolate names",
        color="SIR",
        title="Isolate log2(MIC-values) for different antibiotics",
        range_y=[-11, 12],
        hover_data={
            "Antibiotics": False,
            "Log2(MIC-value)": ":.0f",
            "Scale": True,
            "Pathogen": True,
        },
    )

    # Changes the dot color depending on SIR category
    def change_trace_color(trace):
        if trace.name == "Resistant":
            trace.update(marker_color="tomato")
        elif trace.name == "Intermediate":
            trace.update(marker_color="gold")
        elif trace.name == "Sensitive":
            trace.update(marker_color="limegreen")
        else:
            raise ValueError("Not a valid trace")

    def update_trace(trace):
        """Function that updates traces"""
        change_trace_color(trace)

    # Update dot color
    fig.for_each_trace(update_trace)

    # Modify x-ticks, y-ticks, legend and title
    fig.update_layout(
        xaxis=dict(tickmode="array", tickvals=x_axis, ticktext=antibiotics),
        yaxis=dict(
            tickmode="array",
            tickvals=[i for i in range(-10, 12)],
            ticktext=y_axis_ticktext,
        ),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        title_x=0.5,
    )

    fig.show()


def main():
    # Load files
    chosen_isolates_list = pd.read_csv("Visualisation/Chosen_isolates_list.csv")
    CIB = pd.ExcelFile("Visualisation/CIB_TF-data_AllIsolates_20230302.xlsx")
    matrix_EU = pd.read_excel(CIB, "matrix EU")

    # Rename a long name for plotting purposes
    matrix_EU.rename(
        columns={"Trimethoprim-sulfamethoxazole": "Trimeth-sulf"}, inplace=True
    )

    # Select isolates
    chosen_isolates = extract_chosen_isolates(chosen_isolates_list, matrix_EU)

    # List of antiiotic names
    antibiotics = list(chosen_isolates.columns[3:])

    # Extract all SIRs for an antibiotic.
    chosen_isolates_SIR = extract_SIR(chosen_isolates, antibiotics)

    # Remove the tuples that have None in their SIR data
    filtered_chosen_isolates_SIR = filter_mic_values(chosen_isolates_SIR)

    # Extract the mic-values of each isolate for each antibiotic.
    mic_data = extract_mic_data(filtered_chosen_isolates_SIR, antibiotics)

    # Create dataframe used for plotting
    plot_df = create_plot_df(antibiotics, mic_data)

    plotly_dotplot(plot_df, antibiotics)


if __name__ == "__main__":
    main()
