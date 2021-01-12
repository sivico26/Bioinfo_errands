#! /usr/bin/env python
import pandas as pd
from Bio import SeqIO
from numpy import nan
import argparse

def get_arguments():
    parser = argparse.ArgumentParser(description='''Recalculate the busco score for given subassemblies and 
                                     the Busco full report of the assembly.''',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', type=str, default="./busco_update.tsv",
                        help='Path to put updated scores table (%(default)s)')
    parser.add_argument('-b', '--busco_report', type=str, required=True,
                        help='Path to the full busco report for the assembly')
    parser.add_argument('-i', '--input', type=str, nargs="+",
                        help='Paths to assemblies')
    
    args = parser.parse_args()
    return args

def busco_score(df):
    comp = (df.Status == "Complete").sum()
    frag = (df.Status == "Fragmented").sum()
    miss = (df.Status == "Missing").sum()
    dups = len(set(df[df.Status == "Duplicated"]["# Busco id"]))
    N = sum([comp, dups, frag, miss])
    # C:96.1%[S:91.3%,D:4.8%],F:1.9%,M:2.0%,n:2121
    #print(f"C:{100*(comp+dups)/N :.1f}%[S:{100*comp/N:.1f}%,D:{100*dups/N:.1f}%],F:{100*frag/N:.1f}%,M:{100*miss/N:.1f}%,n:{N}")
    return comp, dups, frag, miss

def busco_update(df, asm):
    busco_orig = list(busco_score(df))
    
    asm_scf = {req.id for req in SeqIO.parse(asm, "fasta")}
    lost_scf = set(df.Contig) - asm_scf
    lost_scf.discard(nan)
    
    if not lost_scf:
        N = sum(busco_orig)        
        print("New Busco score is the same as the original score")
        print("Original Busco Score:")
        print(f"C:{100*(sum(busco_orig[:2]))/N :.1f}%[S:{100*busco_orig[0]/N:.1f}%,D:{100*busco_orig[1]/N:.1f}%],F:{100*busco_orig[2]/N:.1f}%,M:{100*busco_orig[3]/N:.1f}%,n:{N}")
        return tuple(busco_orig)
    else:
        df = df[~(df.Contig.isin(lost_scf))]  ## We filter registers of missing contigs
        dup_counts = df[df.Status == "Duplicated"].groupby(by=["# Busco id"]).count()
        idx = dup_counts[dup_counts.Status == 1].index
        bool_row = (df["# Busco id"].isin(idx)) & (df.Status == "Duplicated")
        
        if bool_row.any():
            df.loc[bool_row,"Status"] = "Complete"     
        
        busco_new = list(busco_score(df))
        busco_new[3] += sum(orig - new for orig, new in zip(busco_orig,busco_new))
        
        N = sum(busco_orig)
        print("Original Busco Score:")
        print(f"C:{100*(sum(busco_orig[:2]))/N :.1f}%[S:{100*busco_orig[0]/N:.1f}%,D:{100*busco_orig[1]/N:.1f}%],F:{100*busco_orig[2]/N:.1f}%,M:{100*busco_orig[3]/N:.1f}%,n:{N}")
        N = sum(busco_new)
        print("New Busco Score:")
        print(f"C:{100*(sum(busco_new[:2]))/N :.1f}%[S:{100*busco_new[0]/N:.1f}%,D:{100*busco_new[1]/N:.1f}%],F:{100*busco_new[2]/N:.1f}%,M:{100*busco_new[3]/N:.1f}%,n:{N}")
        
        return tuple(busco_new)
        
def main():
    args = get_arguments()
    with open(args.output, "w") as out:
        df = pd.read_csv(args.busco_report, sep="\t", skiprows=4, dtype=str)
        for asm in args.input:
            name = asm.split("/")[-1].split(".")[0] + "\t"
            print(f"starting {name}")
            res = name + "\t".join(map(str, busco_update(df, asm))) + "\n"
            out.write(res)
    
if __name__ == "__main__":
    main()  
