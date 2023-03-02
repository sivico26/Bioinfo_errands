#! /usr/bin/env python
from typing import Optional
from Bio import SeqIO
from pathlib import Path
from itertools import accumulate, repeat
from operator import itemgetter, truediv
import matplotlib.pyplot as plt
from matplotlib.colors import TABLEAU_COLORS
import argparse


def get_arguments():
    parser = argparse.ArgumentParser(
        description="Plot ngx contiguity for a set of assemblies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default="./ngx.png",
        help="Path to put contiguity figure (%(default)s)",
    )
    parser.add_argument(
        "-g",
        "--genome_size",
        type=int,
        default=None,
        help="Genome size to calculate NG(X) instead of N(X)",
    )
    parser.add_argument("--ylog", action="store_true", help="Use log scale for y-axis")
    parser.add_argument("--ymin", type=int, default=0, help="Minimum value for y-axis")
    parser.add_argument(
        "--ymax", type=int, default=None, help="Maximum value for y-axis"
    )

    parser.add_argument(
        "-c",
        "--colors",
        type=str,
        nargs="+",
        default=None,
        help="Colors to use in plot in correspondence with assemblies order",
    )
    parser.add_argument(
        "-i", "--input", type=Path, nargs="+", help="Paths to assemblies"
    )

    args = parser.parse_args()
    return args


def NGx(
    output: Path,
    colors: str,
    *assemblies: list[Path],
    genome_size: Optional[int] = None,
    ylog: bool = False,
    ymin: int = 0,
    ymax: Optional[int] = None,
):
    lens = {}
    fraction_len = {}
    for asm in assemblies:
        name = asm.stem
        lens[name] = sorted(
            (len(req) for req in SeqIO.parse(asm, "fasta")), reverse=True
        )

        if not genome_size:
            g = sum(lens[name])
            fraction_len[name] = list(map(truediv, accumulate(lens[name]), repeat(g)))
        else:
            filt = lambda x: x < 1
            fraction_len[name] = list(
                filter(filt, map(truediv, accumulate(lens[name]), repeat(genome_size)))
            )
            lens[name] = lens[name][: len(fraction_len[name])]

    colors = list(TABLEAU_COLORS.values()) if colors is None else colors

    for name, color in zip(lens, colors[: len(lens)]):
        plt.plot(
            [0] + fraction_len[name],
            [lens[name][0]] + lens[name],
            drawstyle="steps-pre",
            color=color,
        )

    if ylog:
        plt.yscale("log")
        plt.legend(list(lens), loc="lower left")
    else:
        plt.legend(list(lens), loc="upper right")

    _, _, yl, ym = plt.axis()
    bottom = ymin if ymin else yl
    top = ymax if ymax else ym
    plt.ylim((bottom, top))

    plt.ylabel("Size (bp)")
    plt.xlabel("Cumulative size fraction")
    plt.title("NG(x)" if genome_size else "N(X)")
    plt.vlines(0.5, 0, top, linestyles="dashed", alpha=0.5)
    plt.savefig(output)


def main():
    args = get_arguments()
    if args.colors != None and len(args.colors) != len(args.input):
        print(
            f"Incompatible colors ({len(args.colors)}) for assemblies ({len(args.input)})"
        )
    else:
        NGx(
            args.output,
            args.colors,
            *args.input,
            genome_size=args.genome_size,
            ylog=args.ylog,
            ymin=args.ymin,
            ymax=args.ymax,
        )


if __name__ == "__main__":
    main()
