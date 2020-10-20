"""
Author: Andrew Harris
Python 3.8
"""
import argparse
from pathlib import Path

import pandas as pd
import plotly
import plotly.express as px


def main():
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
        '-o',
        '--output',
        type=str,
        action='store',
        required=True,
        help='Output plot filename (svg, pdf, png, jpeg).',
    )
    args = parser.parse_args()
    
    # input variables
    INPUT_EXCEL = Path(args.input)
    OUTPUT_PLOT = Path(args.output)

    # read excel file
    excel_data = pd.read_excel(INPUT_EXCEL)
    
    # make plot + output
    fig = px.bar(
        excel_data,
        x="Chromosome",
        y="Length (bp)",
        color="Assembly",
        color_discrete_sequence =["red", "#fdbb84", "#a6bddb"],
        template="simple_white",
        height=500,
        width=1250,
    )
    fig.update_layout(
        barmode='group',
        font=dict(
            family="Arial, monospace",
            size=12,
        ),
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.98,
        )
        )
    fig.show()
    fig.write_image(
        file=OUTPUT_PLOT,
        format=str(INPUT_EXCEL.suffix).strip('.'),
        engine="kaleido",
    )
    return


if __name__ == "__main__":
    main()
