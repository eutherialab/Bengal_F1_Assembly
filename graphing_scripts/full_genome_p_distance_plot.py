"""
Author: Andrew Harris
Python 3.8.3
"""
import argparse
from pathlib import Path

import pandas as pd
import plotly
import plotly.graph_objs as go
import plotly.express as px
from plotly.subplots import make_subplots


def _read_file_to_df(f):
    """Read in file based of suffix"""
    if f.suffix == ".tsv":
        df = pd.read_csv(f, sep='\t')
        return df
    elif f.suffix == '.csv':
        df = pd.read_csv(f, sep=',')
        return df


def make_plot(df, outdir, fname):
    # ---- df melt + organize ----
    df = df.melt(id_vars=["Chromosome", "Start", "Stop"])
    df.columns = ["Chromosome", "Start", "Stop", "Sample", "p_distance"]
    df.sort_values(by=["Chromosome", "Start", "Stop"], inplace=True)

    # ---- collect df info ----
    chromosomes = [i for i in df['Chromosome'].unique()]
    samples = [i for i in df['Sample'].unique()]
    samples.sort(reverse=True)

    # ---- Get max p_distance value to set y-axis range ----
    y_max = float(df["p_distance"].max() * 1.1)  # increase 10% above max
    x_max = df["Stop"].max() * 1.01  # increase 10% above max

    colors = px.colors.qualitative.Prism

    base_graph_height = 150

    fig = make_subplots(
        rows=len(chromosomes),
        cols=1,
        x_title="Position",
        y_title="p-distance",
        row_titles=chromosomes,
        vertical_spacing=0.01,
        row_heights=[base_graph_height] * len(chromosomes),
    )

    for placeholder, sample in enumerate(samples):
        legend_flag = True
        for row, current_chrom in enumerate(chromosomes, start=1):
            filt = (df['Chromosome'] == current_chrom) & (df["Sample"] == sample)
            sample_chromosome_data = df[filt]
            # Make figure
            fig.add_trace(
                go.Scatter(
                    x=sample_chromosome_data['Stop'],
                    y=sample_chromosome_data['p_distance'],
                    mode='lines',
                    legendgroup=str(sample),
                    name=sample,
                    line=dict(
                        color=colors[placeholder],
                        width=float(0.75)
                    ),
                    showlegend=legend_flag,
                ),
                row=row,
                col=1
            )
            legend_flag = False
            continue
        continue
    #  ---- Update figure ----
    fig.update_layout(
        height=base_graph_height*len(chromosomes),
        template="simple_white",
        margin=dict(
            l=60,
            r=10,
            b=60,
            t=0,
            pad=5
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1,
            xanchor="left",
            x=0,
            itemsizing='trace',
            title="",
        ),
        font=dict(
            family="Arial, monospace",
            size=12,
        ),
        annotations=[{
            "font": dict(
                    family="Arial, monospace",
                    size=12,
                    ),
        }]
    )
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig.update_yaxes(range=[0.0, y_max], fixedrange=True)
    fig.update_xaxes(range=[0, x_max])
    html_filename = outdir / f"{fname}.html"
    fig.write_html(str(html_filename))
    return


def main():
    parser = argparse.ArgumentParser(
        description='Creates full genome p-distance plot',
    )
    parser.add_argument(
        '-i', 
        '--input', 
        type=str, 
        action='store', 
        required=True, 
        help='Input p-distance file'
    )
    parser.add_argument(
        '-o', 
        '--output', 
        type=str, 
        action='store', 
        required=True, 
        help='Output location (i.e directory)',
    )
    args = parser.parse_args()
    
    # --- input variables ---
    INPUT_FILE = Path(args.input)
    OUTPUT_DIR = Path(args.output)

    # --- read file + plot ----
    pdist_df = _read_file_to_df(INPUT_FILE)
    make_plot(pdist_df, OUTPUT_DIR, INPUT_FILE.stem)
    return


if __name__ == "__main__":
    main()
