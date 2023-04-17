import pandas as pd
import numpy as np
import plotly.express as px
from data_extraction_functions import (
    extract_chosen_isolates,
    extract_mic_data,
    extract_SIR,
    filter_mic_values,
)


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
    MIC_value_list = []
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
            MIC_value_list.append(2 ** (mic_value))
            on_off_scale.append(scale)
            pathogen_list.append(pathogen)
            if scale is True:
                y_values.append(mic_value_jitter)

            # If off-scale move value to MAX_C or MIN_C
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
            "MIC value": MIC_value_list,
        },
        index=np.arange(len(x_values)),
    )

    return plot_df


def add_rectangles_to_plot(fig, antibiotics: list) -> None:
    names_to_conc = {
        "Benzylpenicillin": ["0.015 - 32.0", "0.008 - 16.0"],
        "Ampicillin": ["0.125 - 64.0", "0.008 - 16.0"],
        "Cefoxitin": ["0.25 - 16.0"],
        "Ceftaroline": ["0.03125 - 16.0"],
        "Ceftobiprole": ["0.03125 - 16.0"],
        "Ceftriaxone": ["0.0078125 - 16.0"],
        "Imipenem": ["0.03125 - 16.0"],
        "Meropenem": ["0.0078125 - 8.0"],
        "Ciprofloxacin": ["0.0625 - 8.0"],
        "Levofloxacin": ["0.03125 - 16.0"],
        "Gentamicin": ["128.0 - 500.0"],
        "Dalbavancin": ["0.015 - 1.0"],
        "Teicoplanin": ["0.03125 - 16.0", "0.03125 - 8.0"],
        "Vancomycin": ["0.125 - 64.0", "0.015 - 8.0"],
        "Erythromycin": ["0.015 - 16.0", "0.004 - 4.0"],
        "Clindamycin": ["0.015 - 16.0", "0.0078125 - 16.0"],
        "Tetracycline": ["0.015 - 32.0", "0.015 - 16.0"],
        "Linezolid": ["0.125 - 16.0", "0.0625 - 8.0"],
        "Daptomycin": ["0.03125 - 16.0"],
        "Rifampicin": ["0.002 - 8.0"],
        "Trimeth-sulf": ["0.0625 - 16.0", "0.0625 - 8.0"],
    }

    for i in range(len(antibiotics)):
        if antibiotics[i] in names_to_conc:
            conc = names_to_conc[antibiotics[i]][0].split("-")
            conc_low = float(conc[0])
            conc_high = float(conc[1])

            box_witdh = 0.4
            fig.add_vrect(
                x0=i - box_witdh,
                x1=i + box_witdh,
                y0=0.99 - (10 - np.log2(conc_low / 4)) / 23,
                y1=1.01 - (10 - np.log2(conc_high / 4)) / 23,
                fillcolor="maroon",
                layer="below",
                line_width=0,
                opacity=0.8,
            )

            if len(names_to_conc[antibiotics[i]]) > 1:
                conc = names_to_conc[antibiotics[i]][1].split("-")
                conc_low = float(conc[0])
                conc_high = float(conc[1])

                box_witdh = 0.3
                fig.add_vrect(
                    x0=i - box_witdh,
                    x1=i + box_witdh,
                    y0=0.99 - (10 - np.log2(conc_low / 4)) / 23,
                    y1=1.01 - (10 - np.log2(conc_high / 4)) / 23,
                    fillcolor="ghostwhite",
                    layer="below",
                    line_width=0,
                    opacity=0.3,
                )

    fig.add_vrect(
        x0=10,
        x1=16,
        y0=0.97,
        y1=1,
        fillcolor="maroon",
        line_width=0,
        opacity=0.8,
        annotation_text="Non-fastidious concentration ranges",
        annotation_position="top",
    )
    fig.add_vrect(
        x0=16.5,
        x1=22.5,
        y0=0.97,
        y1=1,
        fillcolor="ghostwhite",
        line_width=0,
        opacity=0.3,
        annotation_text="Fastidious concentration ranges",
        annotation_position="top",
    )


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
        opacity=0.6,
        title="Isolate MIC-values for different antibiotics",
        range_y=[-11, 12],
        template="plotly_dark",
        hover_data={
            "Antibiotics": False,
            "Log2(MIC-value)": False,
            "Scale": False,
            "Pathogen": True,
            "MIC value": True,
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

    # Update dot color
    fig.for_each_trace(change_trace_color)

    add_rectangles_to_plot(fig, antibiotics)

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
    fig.write_html("first_figure.html", auto_open=True)
    # fig.show()


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
