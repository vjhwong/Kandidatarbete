# Author: Yasmine Sundelin Tjärnström
# Date: Feb/March 2023

# Script for selecting isolates, adapted for student project

import pandas as pd
import json
import numpy as np
import itertools
from read_cib import get_data, rank_system

# Neccessary files:
# read_cib.py
# parameters_settings.py (change parameters here)
# market_prio.json
# ranges.json
# abx_abbr.json
# raw data (input, CIB)

# Isolate selection in multiple steps
# Isolates are selected from a list that is sorted from most to least interesting, according to a point system
# --> most interesting isolates are selected first
# 1. Isolate per species selected
# 2. Isolate selection to fill gaps for each bugdrug combination (with "interesting" isolates)
# 3. Isolate selection to fill to an upper limit (if not already surpassed)


# Features:
# Choose target number of isolates for species (group and subspecies) (selection 1)
# Choose target number of "interesting" isolates (in total) for each bugdrug combination (selection 2)
# Choose which characteristics to look for when filling gaps (selection 2). Choose how many combinations one wants with name "Bugdrug fill isolates x".
# Choose upper fill limit and lower stop limit in isolate selection
# Choose point system of interesting characteristics of isolates
# Other settings such as dataset and software kit version available
# --> See file "parameters_settings.json"


def species_fill(
    chosen_isolates, available_data, tempdata, parameters, isos_req, errors, pat
):

    if parameters["Lower limit"][0]:
        diff = parameters["Lower limit"][1] - len(chosen_isolates)
        if diff < 0:
            chosen_isolates = chosen_isolates[: parameters["Lower limit"][1]]
            return [chosen_isolates, available_data, tempdata, errors]
        if (
            diff < isos_req
        ):  # example: isos_req=5 but we are 2 isolates away from limit --> only choose 2
            isos_req = diff

    chosen_data = tempdata[:isos_req]
    chosen_isolates = pd.concat([chosen_isolates, chosen_data])
    tempdata = tempdata[-tempdata["Isolate"].isin(list(chosen_data["Isolate"]))]
    available_data = available_data[
        -available_data["Isolate"].isin(list(chosen_data["Isolate"]))
    ]

    if len(chosen_data) != isos_req:
        errors = pd.concat(
            [
                errors,
                pd.DataFrame(
                    {
                        "Pathogen": [f"{pat}:"],
                        "Message": f"Not enough isolates in first selection, {len(chosen_data)}/{isos_req} isolates were selected",
                    }
                ),
            ]
        )

    return [chosen_isolates, available_data, tempdata, errors]


def get_bugdrug_fill(parameters):

    # get requirements for bugdrug fill
    bugdrug_fill = list()
    try:
        bugdrug_fill.append(
            [
                parameters[f"Bugdrug fill requirements 1"]["SIR"],
                parameters[f"Bugdrug fill requirements 1"]["scale"],
                parameters[f"Bugdrug fill requirements 1"]["POS"],
            ]
        )
    except:
        print("Need at least one bugdrug fill scenario if bugdrug fill is set to true")

    for i in range(2, 10):

        try:
            bugdrug_fill.append(
                [
                    parameters[f"Bugdrug fill requirements {i}"]["SIR"],
                    parameters[f"Bugdrug fill requirements {i}"]["scale"],
                    parameters[f"Bugdrug fill requirements {i}"]["POS"],
                ]
            )
        except:
            break

    return bugdrug_fill


def bugdrug_fill(
    chosen_isolates, available_data, tempdata, parameters, abx, errors, pat
):

    if not parameters["Bugdrug fill"][0]:
        return [chosen_isolates, available_data, errors]

    else:

        isos_req = parameters["Bugdrug fill"][1]

        scenarios = get_bugdrug_fill(parameters)

        for a in abx:

            remain = isos_req

            if parameters["Lower limit"][0]:
                diff = parameters["Lower limit"][1] - len(chosen_isolates)
                if diff < 0:
                    chosen_isolates = chosen_isolates[: parameters["Lower limit"][1]]
                    return [chosen_isolates, available_data, errors]
                if (
                    diff < remain
                ):  # example: remain=5 but we are 2 isolates away from limit --> only choose 2
                    remain = diff

            # valid data
            isos_valid = [
                row["Isolate"]
                for index, row in tempdata.iterrows()
                if list(row[a].values())[0][0] != 0
            ]
            tempdata_valid = tempdata[tempdata["Isolate"].isin(isos_valid)]

            bugdrug_chosen = chosen_isolates[chosen_isolates["Pathogen"] == pat][[a]]

            for index, row in bugdrug_chosen.iterrows():
                vals = list(row.values[0].values())
                for s in scenarios:
                    # if both US and EU
                    for v in vals:

                        if a == "D-test":
                            POS = True if v[3] == "POS" else False

                            if s[2] != "":
                                if s[2] == POS:  # each can be true or false
                                    remain -= 1
                                    break
                            else:  # req=''
                                if (v[3] == "NEG") | (v[3] == "POS"):  # check if valid
                                    remain -= 1
                                    break

                        # If both EU and US values, only one counts
                        elif (s[0] in str(v[0])) & (
                            s[1] in str(v[3])
                        ):  # SIR and on-scale information
                            remain -= 1
                            break

            for s in scenarios:  # try to find only best scenario first
                for index, row in tempdata_valid.iterrows():
                    # if both US and EU
                    for k, v in row[a].items():
                        if remain < 1:
                            break

                        SIR = v[0]
                        scale = v[3]

                        # v:
                        # [SIR, sign, mic, scale]
                        # [SIR, sign, mic, POS/NEG]

                        # s:
                        # [SIR, scale, POS/NEG]

                        if a == "D-test":

                            POS = True if v[3] == "POS" else False

                            if POS == s[2]:
                                pass  # OK - choose isolates
                            elif s[2] == "":
                                pass  # OK - choose isolates
                            else:
                                continue  # Not OK - move on

                            chosen_data = pd.DataFrame(row.to_frame().T)
                            chosen_isolates = pd.concat([chosen_isolates, chosen_data])
                            tempdata = tempdata[
                                -tempdata["Isolate"].isin(list(chosen_data["Isolate"]))
                            ]
                            available_data = available_data[
                                -available_data["Isolate"].isin(
                                    list(chosen_data["Isolate"])
                                )
                            ]
                            remain -= 1

                        else:
                            if SIR == s[0]:

                                if a == "Gentamicin":
                                    pass  # OK - choose isolates

                                else:
                                    if scale == s[1]:
                                        pass  # OK - choose isolates
                                    else:
                                        continue  # Not OK - move on

                                chosen_data = pd.DataFrame(row.to_frame().T)
                                chosen_isolates = pd.concat(
                                    [chosen_isolates, chosen_data]
                                )
                                tempdata = tempdata[
                                    -tempdata["Isolate"].isin(
                                        list(chosen_data["Isolate"])
                                    )
                                ]
                                available_data = available_data[
                                    -available_data["Isolate"].isin(
                                        list(chosen_data["Isolate"])
                                    )
                                ]
                                remain -= 1

            chosen = isos_req - remain
            if remain > 0:
                if len(isos_valid) < 1:
                    pass
                else:
                    errors = pd.concat(
                        [
                            errors,
                            pd.DataFrame(
                                {
                                    "Pathogen": [f"{a}/{pat}:"],
                                    "Message": f"Not enough interesting isolates, {chosen}/{isos_req} isolates were selected",
                                }
                            ),
                        ]
                    )

    return [chosen_isolates, available_data, errors]


def group_fill(
    chosen_isolates, available_data, parameters, pat_group, errors, subspecies
):

    # fill with other pats from same group if necessary
    if not parameters["Isolates per species"]["Fill group"]:
        return [chosen_isolates, available_data, errors]

    else:

        overall = parameters["Isolates per species"][pat_group]["Overall"]

        # all subspecies within that group
        tempdata = available_data[
            available_data["Pathogen"].isin(
                list(parameters["Isolates per species"][pat_group].keys())[:-2]
            )
        ]

        remain = overall - len(
            chosen_isolates[chosen_isolates["Pathogen"].isin(subspecies)]
        )

        if remain < 0:
            remain = 0

        if parameters["Lower limit"][0]:
            diff = parameters["Lower limit"][1] - len(chosen_isolates)
            if diff < 0:
                chosen_isolates = chosen_isolates[: parameters["Lower limit"][1]]
                return [chosen_isolates, available_data, errors]
            if (
                diff < remain
            ):  # example: remain=5 but we are 2 isolates away from limit --> only choose 2
                remain = diff

        chosen_data = tempdata[:remain]
        chosen_isolates = pd.concat([chosen_isolates, chosen_data])
        available_data = available_data[
            -available_data["Isolate"].isin(list(chosen_data["Isolate"]))
        ]

        # chosen before + chosen now
        chosen_tot = (overall - remain) + len(chosen_data)

        if chosen_tot != overall:
            errors = pd.concat(
                [
                    errors,
                    pd.DataFrame(
                        {
                            "Pathogen": [f"{pat_group}:"],
                            "Message": f"Not enough isolates in group fill, {chosen_tot}/{overall} isolates were selected",
                        }
                    ),
                ]
            )

    return [chosen_isolates, available_data, errors]


def isolate_selection(available_data, parameters, errors, pats_groups, abx):

    chosen_isolates = pd.DataFrame()

    for pat_group in pats_groups:

        subspecies = list(parameters["Isolates per species"][pat_group].keys())[:-2]

        for pat in subspecies:

            isos_req = parameters["Isolates per species"][pat_group][pat]

            tempdata = available_data[available_data["Pathogen"] == pat]

            # Species fill
            [chosen_isolates, available_data, tempdata, errors] = species_fill(
                chosen_isolates,
                available_data,
                tempdata,
                parameters,
                isos_req,
                errors,
                pat,
            )

            # Bugdrug fill
            [chosen_isolates, available_data, errors] = bugdrug_fill(
                chosen_isolates, available_data, tempdata, parameters, abx, errors, pat
            )

        # After going through all subspecies, fill for entire pathogen group
        [chosen_isolates, available_data, errors] = group_fill(
            chosen_isolates, available_data, parameters, pat_group, errors, subspecies
        )

    return [available_data, chosen_isolates, errors]


def upper_fill(available_data, parameters, chosen_isolates):

    # Fill to this limit if not already surpassed

    if (parameters["Upper fill"][0]) & (not parameters["Lower limit"][0]):
        upper_fill = parameters["Upper fill"][1]
        diff = upper_fill - len(chosen_isolates)
        if diff < 0:
            diff = 0
    else:
        return [available_data, chosen_isolates]

    chosen_data = available_data[:diff]
    chosen_isolates = pd.concat([chosen_isolates, chosen_data])
    available_data = available_data[
        -available_data["Isolate"].isin(list(chosen_data["Isolate"]))
    ]

    return [available_data, chosen_isolates]


def iso_sel_setup(available_data, abx, parameters, market_prio):

    # Setup
    errors = pd.DataFrame()
    pat_prio = list(
        dict.fromkeys(
            itertools.chain.from_iterable(
                [market_prio[v].keys() for v, k in market_prio.items()]
            )
        )
    )

    [available_data, chosen_isolates, errors] = isolate_selection(
        available_data, parameters, errors, pat_prio, abx
    )
    [available_data, chosen_isolates] = upper_fill(
        available_data, parameters, chosen_isolates
    )

    return [chosen_isolates, errors]


def main(CIB, parameters, ranges, abx_abbr, market_prio):

    # Setup
    inputdata = pd.ExcelFile(CIB)
    comb_dataset = pd.DataFrame()
    parameters = json.load(open(parameters))
    ranges = json.load(open(ranges))
    abx_abbr = json.load(open(abx_abbr))
    market_prio = json.load(open(market_prio))

    # get data from CIB and relevant antibiotics
    data = {}
    abx = list()
    for d in parameters["Datasets"]:
        try:
            tempdata = pd.read_excel(inputdata, f"matrix {d}")
            data[f"{d}"] = tempdata
        except ValueError:
            print(f"Data region '{d}' not valid, try 'US' or 'EU'")
            continue
        abx += list(data[f"{d}"].columns[3:])
        data_default = data[f"{d}"]
    abx = list(np.unique(abx))

    # Go through data and put into new df with rank
    for j in range(0, len(data_default)):

        # last rows
        if type(data_default.iloc[j, 1]) == float:
            break

        res = {}
        iso = data_default.iloc[j, 0]
        pat = data_default.iloc[j, 1]
        res["Isolate"] = iso
        res["Pathogen"] = pat

        # get fastidious state
        isolates_per_species = dict(
            list(parameters["Isolates per species"].items())[1:]
        )
        break_loop = False
        for k1, v1 in isolates_per_species.items():
            for k2, v2 in v1.items():
                if pat in k2:
                    fast = isolates_per_species[k1]["Fastidious"]
                    break_loop = True
                    break
            if break_loop:
                break
        res["Fastidious"] = fast

        # iterate through every abx(column) for that isolate
        for a in abx:

            # Get data
            final_data = get_data(data, j, a, ranges, abx_abbr, fast, parameters)
            res[a] = final_data

        # Add rank for that isolate
        res = rank_system(res, parameters["Point system"])

        # Add isolate to new dataset
        comb_dataset = pd.concat([comb_dataset, pd.DataFrame([res])])

    # sort isolates by rank
    sorted_dataset = comb_dataset.sort_values("Q-rank", ascending=False)

    # isolate selection
    [chosen_isolates, errors] = iso_sel_setup(
        sorted_dataset, abx, parameters, market_prio
    )

    return [chosen_isolates, sorted_dataset, errors]


if __name__ == "__main__":

    # Input
    CIB = "Isolate Selection Student Project\CIB_TF-data_AllIsolates_20230302.xlsx"
    parameters = "Isolate Selection Student Project\parameters_settings.json"
    ranges = "Isolate Selection Student Project/ranges.json"
    abx_abbr = "Isolate Selection Student Project/abx_abbr.json"
    market_prio = "Isolate Selection Student Project\market_prio.json"

    [chosen_isolates, sorted_dataset, errors] = main(
        CIB, parameters, ranges, abx_abbr, market_prio
    )

    # Output
    chosen_isolates["Isolate"].to_csv("Chosen_isolates_list.csv", index=False)
    np.savetxt("Errors.txt", errors.to_numpy(), fmt="%s")
