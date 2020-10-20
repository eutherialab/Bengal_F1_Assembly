"""
Author: Andrew Harris
Python Version 3.8
"""
import logging
from pathlib import Path
import statistics

import plotly
import plotly.express as px
import pandas as pd

Path("../logs/").mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename='../logs/side_by_side_plot.log',
                    level=0, filemode='w', format='')


def replace_to(value):
    new_value = value.replace("2", "-to-")
    return new_value


def make_median_list(df):
    median_list = []
    for rl, c in zip(df["Read Length"], df["Count"]):
        for _ in range(c):
            median_list.append(int(rl))
            continue

    median_list.sort()
    return median_list


def build_median_graph_vlines(full_data_df):
    lines = []
    crosses = [c for c in full_data_df["Cross"].unique()]
    subtypes = [
        # "correctly-sorted",
        "mom-to-dad",
        "dad-to-mom",
        "mom-to-unknown",
        "dad-to-unknown",
        "unknown-to-mom",
        "unknown-to-dad",
    ]
    logging.info(f"Crosses: {crosses}")

    cross1 = crosses[0]
    cross2 = crosses[1]

    # xref-labels go from bottom to top (1-4), so need to assign x_axis backwards
    x_axis = len(crosses) * len(subtypes)

    for subtype in subtypes:
        cross1_data = full_data_df[(full_data_df["Cross"] == cross1) & (
            (full_data_df["Subtype"] == subtype))]
        sorted_cross1_data = make_median_list(cross1_data)
        cross1_median = statistics.median(sorted_cross1_data)

        cross2_data = full_data_df[(full_data_df["Cross"] == cross2) & (
            (full_data_df["Subtype"] == subtype))]
        sorted_cross2_data = make_median_list(cross2_data)
        cross2_median = statistics.median(sorted_cross2_data)

        cross1_logging_info = f"Cross: {cross1}, Subtype: {subtype}, Median: {cross1_median}, DF_len: {len(cross1_data)}"
        cross2_logging_info = f"Cross: {cross2}, Subtype: {subtype}, Median: {cross2_median}, Df_len: {len(cross2_data)}"
        logging.info(cross1_logging_info)
        logging.info(cross2_logging_info)

        lines.append(
            dict(
                type="line",
                yref=f'y{x_axis - 1}',
                y0=0,
                y1=cross1_data["Count"].max(),
                xref=f'x{x_axis - 1}',
                x0=cross1_median,
                x1=cross1_median,
                opacity=0.8,
                line=dict(
                    color="black",
                    width=1,
                )
            )
        )
        lines.append(
            dict(
                type="line",
                yref=f'y{x_axis}',
                y0=0,
                y1=cross2_data["Count"].max(),
                xref=f'x{x_axis}',
                x0=cross2_median,
                x1=cross2_median,
                opacity=0.8,
                line=dict(
                    color="black",
                    width=1,
                )
            )
        )
        x_axis -= 2
    return lines


def make_side_by_side_plot():
    """
    This script takes in read length distribution data from two crosses
    and plots them side-by-side, faceted by incorrectly sorted subtypes
    """
    # ---- Input variables ----
    cross_1 = "LilBubxPbe53"
    file_1 = "../misc_data/LilBubxPbe53_read_length_dist_data.tsv"

    cross_2 = "Fca508xPbe14"
    file_2 = "../misc_data/Fca508xPbe14_read_length_dist_data.tsv"

    # Make output directory for plot if not already made
    Path("../plots/").mkdir(parents=True, exist_ok=True)

    # Load in files into df
    df_1 = pd.read_csv(file_1, sep="\t")
    df_2 = pd.read_csv(file_2, sep="\t")

    # Add cross names to dataframes
    df_1["Cross"] = [str(cross_1)] * len(df_1)
    df_2["Cross"] = [str(cross_2)] * len(df_2)

    # Concatenate the dataframes
    concat_df = pd.concat([df_1, df_2])

    # replace '2' with '-to-'
    concat_df["Subtype"] = concat_df["Subtype"].apply(replace_to)

    # Rename columns
    concat_df.columns = ['Read Length', 'Count', 'Subtype', 'Cross']

    # Create median line shapes -- NOT USED
    # fig_median_lines = build_median_graph_vlines(concat_df)

    # Create plot and update settings
    fig = px.scatter(
        concat_df,
        x="Read Length",
        y="Count",
        color="Subtype",
        size_max=2,
        facet_col="Cross",
        facet_row="Subtype",
        facet_row_spacing=0.03,
        template="simple_white",
        height=800,
        width=1000,
        render_mode='svg',
        category_orders={
            "Subtype": [
                "mom-to-dad",
                "dad-to-mom",
                "mom-to-unknown",
                "dad-to-unknown",
                "unknown-to-mom",
                "unknown-to-dad",
            ],
        },
    )
    fig.update_layout(
        legend=dict(
            orientation="h",
            yanchor="top",
            y=-0.1,
            xanchor="center",
            x=0.5,
            title="",
        ),
        font=dict(
            family="Arial, monospace",
            size=12,
        ),
        margin=dict(
            l=5,
            r=5,
        ),
        # shapes=fig_median_lines,
    )
    fig.update_traces(marker=dict(size=4))
    fig.update_xaxes(
        tickmode='linear',
        tick0=0,
        dtick=20000,
        range=[0, 160000]
    )
    fig.update_yaxes(matches=None)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    # fig.show()
    fig.write_image(
        f"../plots/{cross_1}_{cross_2}_Read_Length_Distribution.svg",
        format='svg',
        engine="kaleido",
    )
    return


if __name__ == "__main__":
    make_side_by_side_plot()
