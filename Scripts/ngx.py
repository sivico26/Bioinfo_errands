#! /usr/bin/env python
from Bio import SeqIO
from itertools import accumulate
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def get_arguments():
    parser = argparse.ArgumentParser(description='Plot ngx contiguity for a set of assemblies',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', type=str, default="./ngx.png",
                        help='Path to put contiguity figure (%(default)s)')
    parser.add_argument('-g', '--genome_size', type=int, default=None,
                        help='Genome size to calculate NG(X) instead of N(X)')
    parser.add_argument('--ylog', action ="store_true",
                       help='Use log scale for y-axis')
    parser.add_argument('-c', '--colors', type=str, nargs="+", default=None,
                        help='Colors to use in plot in correspondence with assemblies order')
    parser.add_argument('-i', '--input', type=str, nargs="+",
                        help='Paths to assemblies')
    
    args = parser.parse_args()
    return args

def NGx(output, colors,*assemblies, genome_size = None, ylog = False):
    lens = {}
    fraction_len = {}
    for asm in assemblies:
        name = asm.split("/")[-1].split(".")[0]
        lens[name] = [len(req) for req in SeqIO.parse(asm, "fasta")]
        lens[name].sort(reverse = True)
        
        if not genome_size:
            g = sum(lens[name])
            fraction_len[name] = list(map(lambda x: x/g, accumulate(lens[name]))) + [0]
        else:
            filt = lambda x: x < 1
            norm = lambda x: x/genome_size
            fraction_len[name] = list(filter(filt, map(norm, it.accumulate(lens[name])))) # + [0]    
            lens[name] = lens[name][:len(fraction_len[name])]
            fraction_len[name].append(0) 
    
    #colors = colors[0].split(",")
    longest = 0
    if colors:
        for name, color in zip(lens.keys(), colors):
            lens[name].append(lens[name][0])
            longest = lens[name][0] if lens[name][0] > longest else longest
            sns.lineplot(fraction_len[name], lens[name], drawstyle='steps-pre', color = color)
    else:
        for name in lens.keys():
            lens[name].append(lens[name][0])
            longest = lens[name][0] if lens[name][0] > longest else longest
            sns.lineplot(fraction_len[name], lens[name], drawstyle='steps-pre')

    
    if ylog: 
        plt.yscale("log")
        plt.legend(list(lens.keys()), loc="lower left")
    else:
        plt.legend(list(lens.keys()), loc="upper right")
    
    plt.ylabel("Size (bp)")
    plt.xlabel("Cumulative size fraction")
    plt.title("NG(x)" if genome_size else "N(X)")
    plt.vlines(0.5, 0, longest, linestyles="dashed", alpha = 0.5)
    plt.savefig(output)
    plt.show()

def main():
    args = get_arguments()
    if args.colors != None and len(args.colors) != len(args.input):
        print(f"Incompatible colors ({len(args.colors)}) for assemblies ({len(args.input)})")
    else:
        NGx(args.output, args.colors, *args.input, genome_size = args.genome_size, ylog= args.ylog)

if __name__ == "__main__":
    main()
