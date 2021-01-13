# Bioinformatics errands
A diversity of scripts to solve specific bioinformatics problems.

Below there is a brief description of each script.

## Scripts

### Plotter

`Plotter` was devised to visualize how a given characteristic of a genome vary along the sequences (*e.g* GC content, SNP density, etc.). 

The main input is a two column tabular file where the first column is the genomic position and the second is the value of the characteristic of interest (*e.g.* GC content). 

If multiple sequences are used (several chromosomes, scaffolds, etc.), an additional third column is expected containing the name of the sequence. Plotter concatenates the genomic positions of the different sequences along the x-axis (*i.e* the start of the second chromosome goes just after the end of the first and so on) and adds a vertical with the name of the chromosome to highlight the threshold between sequences.

The script assumes that the tabular file is sorted in both sequence name and genomic positions.

In essence, it just does an scatterplot from a tabular file and performs some formatting. 

```
usage: plotter.py [-h] [-i INPUT] [-o OUTPUT] [-s SIZE SIZE] [--ylab YLAB] [--print-input-examples] [-m]

Plot a basic scatterplot from a tabular file.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to input tabular file. (default: None)
  -o OUTPUT, --output OUTPUT
                        Path to put output image (.png). Default is to show the figure. (default: None)
  -s SIZE SIZE, --size SIZE SIZE
                        Specify the size of the figure X and Y (e.g. 10 5). (default: (13, 8))
  --ylab YLAB           Specify the text to put in the y-axis (e.g. '%GC') (default: None)
  --print-input-examples
                        Print examples of input files (default: None)
  -m, --multi           Specify if the input has multiple sequences. In which case we expect an extra column indicating the name of the sequence. (default: False)
```

### NGx

`NGx` was devised, as its name suggests, to make NGx plots to help to compare different assemblies of the same genome on a contiguity basis. The rationale behind is that metrics such as N50 do not capture the contiguity well enough and is better to visualize the entire distribution or cumulative curve of sequence sizes.

The expected input are the fasta files of the different assemblies to compare.

```
usage: ngx.py [-h] [-o OUTPUT] [-g GENOME_SIZE] [--ylog] [-c COLORS [COLORS ...]] [-i INPUT [INPUT ...]]

Plot ngx contiguity for a set of assemblies

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Path to put contiguity figure (./ngx.png)
  -g GENOME_SIZE, --genome_size GENOME_SIZE
                        Genome size to calculate NG(X) instead of N(X) (default: None)
  --ylog                Use log scale for y-axis (default: False)
  -c COLORS [COLORS ...], --colors COLORS [COLORS ...]
                        Colors to use in plot in correspondence with assemblies order (default: None)
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Paths to assemblies (default: None)
```

### Busco_update

The purpose of `Busco_update` is to update the BUSCO score of an assembly after filtering some of its contigs/scaffolds (see `filter_asm` below).

Its inputs are the full annotation table (generated by busco) and the filtered assembly to update.

```
usage: busco_update.py [-h] [-o OUTPUT] -b BUSCO_REPORT [-i INPUT [INPUT ...]]

Recalculate the busco score for given subassemblies and the Busco full report of the assembly.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Path to put updated scores table (./busco_update.tsv)
  -b BUSCO_REPORT, --busco_report BUSCO_REPORT
                        Path to the full busco report for the assembly (default: None)
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Paths to assemblies (default: None)
```

### Filter_asm

`Filter_asm`, as its name suggests, filters short contigs/scaffolds from an assembly with one of several criteria: 

- **Number:** Keep the N longest sequences from the assembly and filter the remaining.
- **Size:** Filter the sequences shorter than a given size from the assembly.
- **Genome size:** Filter the shortest sequences from an assembly such as the assembly is the closest to a given genome size (slightly larger).

```
usage: filter_asm.py [-h] [-i INPUT] [-o OUTPUT] (-g GENOME_SIZE | -s SIZE | -n NUMBER)

Filter shortest sequences of an assembly

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Paths to assembly (default: None)
  -o OUTPUT, --output OUTPUT
                        Path to filtered assembly (./asm.fasta)
  -g GENOME_SIZE, --genome_size GENOME_SIZE
                        Filter the shortest sequences after genome size is reached (default: None)
  -s SIZE, --size SIZE  Filter sequences smaller than size (default: None)
  -n NUMBER, --number NUMBER
                        Keep the N longest sequences (default: None)
```

### Extract_cds

`Extract_cds` was devised for a quite specific phylogenomics task: given a genbank file, extract the cds sequences and create a fasta file for each cds. If the file is a multigenbank (*e.g.*  several taxa) it outputs a folder for each cds which includes the sequences for each taxa.

```
usage: extract_cds.py [-h] [-i KEY] [-m] input output

Extract CDS from a genbank to output a fasta

positional arguments:
  input                 Path to input genbank file
  output                Path to put file/folder output

optional arguments:
  -h, --help            show this help message and exit
  -i KEY, --ignore KEY  if 'key' matches a CDS name it won't be included in the output (default: None)
  -m, --multi           Specify if the input file is a multigenbank, in which case the CDS of each entry would be extracted in a different fasta file in an output
                        directory in the specified output path (default: False)
```

### Mkblast2

`Mkblast2` uses Biopython blast utility to make blast searches (over the web) of a set of sequences given and parse the hits and summarizes it in tabular format. It can also just generate the summary with the blast hits if they are already produced.

```
usage: mkblast2.py [-h] [-f STR] [-l INT] [-u STR] [-p] [-d STR] [-m INT] input output

Make blasts for sequences from a file and summarise them

positional arguments:
  input                 Path to file input
  output                Path to put output

optional arguments:
  -h, --help            show this help message and exit
  -f STR, --format STR  Input file format [fasta]
  -l INT, --longest INT
                        Whether to take just the N longest sequences in file [None]
  -u STR, --utility STR
                        Blast utility to use [blastn]
  -p, --parse           Do not make the blast searches, only generate summary [False]
  -d STR, --database STR
                        Database to search in [nt]
  -m INT, --matches INT
                        Number of hits to report in summary per sequence [5]

Developed by Simón Villanueva Corrales @sivico26
```

### Clust_filt

`Clust_filt` was designed to take the output of [Philoligo](https://github.com/itsmeludo/PhylOligo) and separate the original assembly to keep the desired clusters found. 

Its inputs are the .dat file generated by `Philoligo` (with cluster information for the sequences) and the assembly that is going to be filtered.

```
usage: clust_filt.py [-h] [-i INPUT] [-d DAT] [-c CLUSTERS [CLUSTERS ...]] [-o OUTPUT] [-f STR]

Filter a genome based on a set of given clusters

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to genome input (default: None)
  -d DAT, --dat DAT     Path to clusters file(.dat) (default: None)
  -c CLUSTERS [CLUSTERS ...], --clusters CLUSTERS [CLUSTERS ...]
                        Spaced list of clusters to be included (default: None)
  -o OUTPUT, --output OUTPUT
                        Path to put output genome [genome_filt.fasta]
  -f STR, --format STR  Input file format [fasta]

Developed by Simón Villanueva Corrales @sivico26
```



