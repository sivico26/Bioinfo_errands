#!/usr/bin/env python
import matplotlib.pyplot as plt
import argparse
import sys

def get_arguments():
    parser = argparse.ArgumentParser(description='Plot a basic scatterplot from a tabular file.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input', type=str, help='Path to input tabular file.')
    parser.add_argument('-o', '--output', type=str, default=None, help='Path to put output image (.png). Default is to show the figure.')
    parser.add_argument('-s', '--size', type=int, nargs=2, default = (13,8), help = "Specify the size of the figure X and Y (e.g. 10 5).")
    parser.add_argument('--ylab', type=str, default = None, help = "Specify the text to put in the y-axis (e.g. '%%GC')")
    parser.add_argument('--print-input-examples', action ='store_true', default=None, help="Print examples of input files")
    parser.add_argument('-m', '--multi', action ='store_true', default = False, help = "Specify if the input has multiple sequences. In which case we expect an extra column indicating the name of the sequence.")

    args = parser.parse_args()
    return args

def parse_tab_input(input_file, multi = True):
    with open(input_file,"r") as tab:
        if multi:
            pos, values, seqs = zip(*(line.split() for line in tab))
            seqs = iter(seqs)
        else:
            pos, values = zip(*(line.split() for line in tab))
        pos = map(int, pos)
        values = map(lambda x: float(x.replace(",",".")), values)
        #values = (float(val.replace(",",".")) for val in values)
        if multi:
            return pos, values, seqs
        else:
            return pos, values
        
def print_examples():
    msg = """We expect a tabular input with 2 columns: Position and value. For instance:
    
        100\t0.38
        200\t0.41
        300\t0.33
        400\t0.36
        
If the 'multi' flag is used we expect an additional column indicating the name of the sequences:
    
        100\t0.38\tchr1
        200\t0.41\tchr1
        300\t0.33\tchr1
        100\t0.45\tchr2
        200\t0.47\tchr2
        """
    print(msg)

def main():
    args = get_arguments()
    multi = args.multi
    output = args.output
    ylab = args.ylab
    size = args.size
    if args.print_input_examples:
        print_examples()
        sys.exit()
    elif not args.input:
        sys.exit("An input to work with is required.")

    if not multi:
        a,b = parse_tab_input(args.input, multi)
        a,b = list(a), list(b)
        mean_val = sum(b)/len(b)
        fig, ax = plt.subplots(figsize=size)
        plt.scatter(x=a, y=b, alpha = 0.5)
        plt.hlines(mean_val, a[0], a[-1],colors="red", linestyles="dashed")
        plt.xlabel("Position (bp)")
        if ylab: plt.ylabel(ylab)
        if output:
            plt.savefig(output, format= "png")
        else:
            plt.show()
    else:
        #Leemos el dataset
        a,b,c = parse_tab_input(args.input, multi)
        a,b,c = list(a), list(b), list(c)
        
        # identificamos las posiciones donde cambian las secuencias
        names = list(set(c))
        if len(names) == 0:
            sys.exit("Names of sequences missing. Are you sure it is a multi-sequence input?")
        names.sort(key=lambda x:c.index(x))
        idx = [c.index(name) for name in names]

        # Sacamos slices para cada secuencia
        s = [slice(j, idx[i+1]) for i,j in enumerate(idx[:-1])] + [slice(idx[-1], None)]

        # Modificamos las posiciones del dataset para concatenar cada secuencia una tras otra
        
        #window_size = a[1] - a[0]
        mod = 0
        new_pos = a[s[0]]
        for i,sli in enumerate(s[1:],1): 
            mod += a[idx[i] - 1]
            new_pos = new_pos + list(map(lambda x: x + mod , a[sli]))

        # Graficamos
        fig, ax = plt.subplots(figsize=size)
        plt.scatter(x=new_pos, y=b, alpha = 0.5)
        # Agregamos los promedios de %GC para cada secuencia
        for sli in s: 
            
            dom = new_pos[sli]
            rang = b[sli]
            mean_val = sum(rang)/len(rang)
            
            plt.hlines(mean_val, dom[0], dom[-1],colors="red", linestyles="dashed")
        # Agregamos los límites donde empieza cada secuancia
        ymin, ymax = min(b), max(b)
        for i in idx:
            plt.vlines(new_pos[i], ymin, ymax, colors="green", label=c[i], linestyles="dashdot")
            plt.annotate(c[i], xy=(new_pos[i], ymax))
        plt.xlabel("Concatenated position (bp)")
        if ylab: plt.ylabel(ylab)
        if output:
            plt.savefig(output, format= "png")
        else:
            plt.show()
        
if __name__ == '__main__':
        sys.exit(main())
