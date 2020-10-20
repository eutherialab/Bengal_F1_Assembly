"""
Author: Andrew Harris
Python 3.8.3
Description: This script plots read count data in a bar plot
Reference: Figure 2 and Supplementary Figure 9
"""
import argparse

import pandas as pd
import plotly
import plotly.express as px


def make_graph():
    # --- Argparse Variables ---
    parser = argparse.ArgumentParser(
        description='Plot read count data for multiple replacement crosses',
    )
    parser.add_argument(
        '-i',
        '--input',
        type=str,
        action='store',
        required=True,
        help='Pathway to input excel file',
    )
    parser.add_argument(
        '-s',
        '--sheet',
        type=str,
        action='store',
        required=True,
        help='Indicate the sheet from the excel file to use.',
    )
    args = parser.parse_args()
    
    input_file = args.input
    sheet = args.sheet
    
    # input_file = "G:\\My Drive\\Projects\\TrioCanu_NonBio_Parent_Project\\excel_files\\Coverage_Info\\Single_Bio_Parent_Cross_just_read_count-draft4.xlsx"
    # sheet = "All_Samples"

    # input_file = "G:\\My Drive\\Projects\\TrioCanu_NonBio_Parent_Project\\excel_files\\Coverage_Info\\Single_Bio_Parent_Cross_just_read_count-draft5.xlsx"
    # sheet = "Figure1C"

    # Read input into a dataframe and sort values
    read_input = pd.read_excel(input_file, sheet_name=sheet)
    read_input.sort_values(by="Read Counts", inplace=True)

    fig = px.bar(
        read_input,
        x="Cross",
        y="Read Counts",
        text="Read Counts",
        color="Parent",
        color_discrete_map={"Maternal Haplotype": "red",
                            "Paternal Haplotype": "blue"},
        height=500,
        width=1250,
        template="simple_white",
        # category_orders={
        #     "Cross": [
        #         "Fca508*xPbe38",
        #         "Fca508*xPbe14",
        #         "Fca508*xPbe53*",
        #         "LilBubxPbe53*",
        #         "PammixPbe53*",
        #         "Fca508*xPbe6",
        #     ],
        # },
        category_orders={
            "Cross": [
                "Fca508* x Pbe53*",
                "CatIIIxPbe53*",
                "LilBub x Pbe53*",
                "Pammi x Pbe53*",
                "Fca508* x Pbe14",
                "Fca508* x Pbe38",
                "Fca508* x Pbe6",
                "Pammi x Pvi12",
                "Fca508* x Pja5",
                "Fca508* x Pja25",
            ],
        },
        # category_orders={
        #     "Cross": [
        #         "Fca508* x Pbe14",
        #         "Fca508* x Pbe53*",
        #         "LilBub x Pbe53*",
        #     ],
        # },
    )
    fig.update_layout(
        barmode='group',
        font=dict(
            family="Arial, monospace",
            size=12,
        ),
        xaxis=dict(
            tickangle=0,
            title_standoff=25,
        ),
        yaxis=dict(
            range=[0, 8000000]
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
            itemsizing='trace',
        ),
        legend_title_text='',
    )
    fig.update_traces(texttemplate='%{text:.2s}', textposition='outside')
    fig.show()
    fig.write_image(file="../plots/Sup_Fig9_500x1250px.svg",
                    format="svg", engine="kaleido")
    return


if __name__ == "__main__":
    make_graph()
