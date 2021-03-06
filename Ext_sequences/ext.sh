#!/bin/bash

# We need a fasta with a query or several --query
# A bool that tells if there is a multi query? Grep might make this unnecesary
# A reference sequence --ref
# An option to indicate splicing or not queries
# UTRs
# Upstream input (bp) -u
# Downstream input (bp) -l
# Normal output would be only exons
# Whole flag to indicate to colapse the alignment

# ARG_OPTIONAL_SINGLE([query],[q],[Path to sequence to align to reference (.fasta).])
# ARG_OPTIONAL_SINGLE([alignment],[a],[Path to alignment (.sam). If provided --query and --out-aln would be ignored.])
# ARG_OPTIONAL_SINGLE([reference],[r],[Path to reference to align to (.fasta).])
# ARG_OPTIONAL_SINGLE([out-aln],[o],[Path to save resulting alignment (.sam).])
# ARG_OPTIONAL_SINGLE([sequence],[s],[Path to save extracted sequences (.fasta).])
# ARG_OPTIONAL_BOOLEAN([exons],[e],[Split the intervals found in the alignment and extract those sequences. Not compatible with --utrs. [false].],[off])
# ARG_OPTIONAL_BOOLEAN([no-splice],[n],[Turn off the splice awareness of the alignment [false].],[off])
# ARG_OPTIONAL_BOOLEAN([gene-exclusion],[g],[Exclude the main alignment from the extracted sequences. Use with --utrs to extract only UTRs [false].],[off])
# ARG_OPTIONAL_BOOLEAN([utrs],[t],[Include upstream and downstream bases of the alignment in the extracted sequences [false].],[off])
# ARG_OPTIONAL_SINGLE([upbases],[u],[Number of upstream bases to extract [200].],[200])
# ARG_OPTIONAL_SINGLE([downbases],[d],[Number of downstream bases to extract [200].],[200])
# ARG_HELP([Usage:  Ext.sh -r reference -s sequence [ -q query & -o out_alignment | -a alignment ] [ --no-splice ] [ --exons ] [ --gene-exclusion ] [--utrs ] [-u upbases] [-d downbases]
#Extract sequences from a reference given an alignment or another sequence that align with it.
#Must have magicblast, bioawk, and bedtools available on the PATH])
# ARGBASH_GO()
# needed because of Argbash --> m4_ignore([
### START OF CODE GENERATED BY Argbash v2.9.0 one line above ###
# Argbash is a bash code generator used to get arguments parsing right.
# Argbash is FREE SOFTWARE, see https://argbash.io for more info
# Generated online by https://argbash.io/generate


die()
{
	local _ret="${2:-1}"
	test "${_PRINT_HELP:-no}" = yes && print_help >&2
	echo "$1" >&2
	exit "${_ret}"
}


begins_with_short_option()
{
	local first_option all_short_options='qarosengtudh'
	first_option="${1:0:1}"
	test "$all_short_options" = "${all_short_options/$first_option/}" && return 1 || return 0
}

# THE DEFAULTS INITIALIZATION - OPTIONALS
_arg_query=
_arg_alignment=
_arg_reference=
_arg_out_aln=
_arg_sequence=
_arg_exons="off"
_arg_no_splice="off"
_arg_gene_exclusion="off"
_arg_utrs="off"
_arg_upbases="200"
_arg_downbases="200"


print_help()
{
#	printf '%s\n' "Usage:  Ext.sh -r reference -s sequence [ -q query & -o out_alignment | -a alignment ] [ --no-splice ] [ --exons ] [ --gene-exclusion ] [--utrs ] [-u upbases] [-d downbases]
	printf '%s\n\n' "Extract sequences from a reference given an alignment or another sequence that alings to it."
	printf 'Usage: %s [-q|--query <arg>] [-a|--alignment <arg>] [-r|--reference <arg>] [-o|--out-aln <arg>] [-s|--sequence <arg>] [-e|--(no-)exons] [-n|--(no-)no-splice] [-g|--(no-)gene-exclusion] [-t|--(no-)utrs] [-u|--upbases <arg>] [-d|--downbases <arg>] [-h|--help]\n\n' "$0"
	printf '\t%s\n' "-q, --query: Path to sequence to align to reference (.fasta). (no default)"
	printf '\t%s\n' "-a, --alignment: Path to alignment (.sam). If provided --query and --out-aln would be ignored. (no default)"
	printf '\t%s\n' "-r, --reference: Path to reference to align to (.fasta). (no default)"
	printf '\t%s\n' "-o, --out-aln: Path to save resulting alignment (.sam). (no default)"
	printf '\t%s\n' "-s, --sequence: Path to save extracted sequences (.fasta). (no default)"
	printf '\t%s\n' "-e, --exons, --no-exons: Split the intervals found in the alignment and extract those sequences. Not compatible with --utrs. [false]. (off by default)"
	printf '\t%s\n' "-n, --no-splice, --no-no-splice: Turn off the splice awareness of the alignment [false]. (off by default)"
	printf '\t%s\n' "-g, --gene-exclusion, --no-gene-exclusion: Exclude the main alignment from the extracted sequences. Use with --utrs to extract only UTRs [false]. (off by default)"
	printf '\t%s\n' "-t, --utrs, --no-utrs: Include upstream and downstream bases of the alignment in the extracted sequences [false]. (off by default)"
	printf '\t%s\n' "-u, --upbases: Number of upstream bases to extract [200]. (default: '200')"
	printf '\t%s\n' "-d, --downbases: Number of downstream bases to extract [200]. (default: '200')"
	printf '\t%s\n\n' "-h, --help: Prints help"
	printf '%s\n' "Must have magicblast, bioawk, and bedtools available on the PATH."
}


parse_commandline()
{
	while test $# -gt 0
	do
		_key="$1"
		case "$_key" in
			-q|--query)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_query="$2"
				shift
				;;
			--query=*)
				_arg_query="${_key##--query=}"
				;;
			-q*)
				_arg_query="${_key##-q}"
				;;
			-a|--alignment)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_alignment="$2"
				shift
				;;
			--alignment=*)
				_arg_alignment="${_key##--alignment=}"
				;;
			-a*)
				_arg_alignment="${_key##-a}"
				;;
			-r|--reference)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_reference="$2"
				shift
				;;
			--reference=*)
				_arg_reference="${_key##--reference=}"
				;;
			-r*)
				_arg_reference="${_key##-r}"
				;;
			-o|--out-aln)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_out_aln="$2"
				shift
				;;
			--out-aln=*)
				_arg_out_aln="${_key##--out-aln=}"
				;;
			-o*)
				_arg_out_aln="${_key##-o}"
				;;
			-s|--sequence)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_sequence="$2"
				shift
				;;
			--sequence=*)
				_arg_sequence="${_key##--sequence=}"
				;;
			-s*)
				_arg_sequence="${_key##-s}"
				;;
			-e|--no-exons|--exons)
				_arg_exons="on"
				test "${1:0:5}" = "--no-" && _arg_exons="off"
				;;
			-e*)
				_arg_exons="on"
				_next="${_key##-e}"
				if test -n "$_next" -a "$_next" != "$_key"
				then
					{ begins_with_short_option "$_next" && shift && set -- "-e" "-${_next}" "$@"; } || die "The short option '$_key' can't be decomposed to ${_key:0:2} and -${_key:2}, because ${_key:0:2} doesn't accept value and '-${_key:2:1}' doesn't correspond to a short option."
				fi
				;;
			-n|--no-no-splice|--no-splice)
				_arg_no_splice="on"
				test "${1:0:5}" = "--no-" && _arg_no_splice="off"
				;;
			-n*)
				_arg_no_splice="on"
				_next="${_key##-n}"
				if test -n "$_next" -a "$_next" != "$_key"
				then
					{ begins_with_short_option "$_next" && shift && set -- "-n" "-${_next}" "$@"; } || die "The short option '$_key' can't be decomposed to ${_key:0:2} and -${_key:2}, because ${_key:0:2} doesn't accept value and '-${_key:2:1}' doesn't correspond to a short option."
				fi
				;;
			-g|--no-gene-exclusion|--gene-exclusion)
				_arg_gene_exclusion="on"
				test "${1:0:5}" = "--no-" && _arg_gene_exclusion="off"
				;;
			-g*)
				_arg_gene_exclusion="on"
				_next="${_key##-g}"
				if test -n "$_next" -a "$_next" != "$_key"
				then
					{ begins_with_short_option "$_next" && shift && set -- "-g" "-${_next}" "$@"; } || die "The short option '$_key' can't be decomposed to ${_key:0:2} and -${_key:2}, because ${_key:0:2} doesn't accept value and '-${_key:2:1}' doesn't correspond to a short option."
				fi
				;;
			-t|--no-utrs|--utrs)
				_arg_utrs="on"
				test "${1:0:5}" = "--no-" && _arg_utrs="off"
				;;
			-t*)
				_arg_utrs="on"
				_next="${_key##-t}"
				if test -n "$_next" -a "$_next" != "$_key"
				then
					{ begins_with_short_option "$_next" && shift && set -- "-t" "-${_next}" "$@"; } || die "The short option '$_key' can't be decomposed to ${_key:0:2} and -${_key:2}, because ${_key:0:2} doesn't accept value and '-${_key:2:1}' doesn't correspond to a short option."
				fi
				;;
			-u|--upbases)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_upbases="$2"
				shift
				;;
			--upbases=*)
				_arg_upbases="${_key##--upbases=}"
				;;
			-u*)
				_arg_upbases="${_key##-u}"
				;;
			-d|--downbases)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_downbases="$2"
				shift
				;;
			--downbases=*)
				_arg_downbases="${_key##--downbases=}"
				;;
			-d*)
				_arg_downbases="${_key##-d}"
				;;
			-h|--help)
				print_help
				exit 0
				;;
			-h*)
				print_help
				exit 0
				;;
			*)
				_PRINT_HELP=yes die "FATAL ERROR: Got an unexpected argument '$1'" 1
				;;
		esac
		shift
	done
}

parse_commandline "$@"

# Argbash ends here.

if [[ $# -eq 0 ]] ; then
	print_help
	exit 0
fi

QUERY=$_arg_query
if [ -z $_arg_alignment ] ; then ALN=false ; else ALN=$_arg_alignment ; fi
#ALN=$_arg_alignment
REF=$_arg_reference
OUTALN=$_arg_out_aln
SEQ=$_arg_sequence
EXONS=$_arg_exons
if [ $_arg_no_splice = "on" ]; then SPLICE=false ; else SPLICE=true; fi
#SPLICE=$_arg_no_splice
UTRS=$_arg_utrs
UPBASES=$_arg_upbases
DOWNBASES=$_arg_downbases
if [ $_arg_gene_exclusion = "on" ]; then GENEXC=flank ; else GENEXC=slop; fi

#GENEXC=$_arg_gene_exclusion

# REF # Path a la referencia a la cual se alineará (.fasta)
# ALN # archivo con alineamiento (.sam/.bam)
# OUTALN # Path donde colocar el alineamiento (.sam)
# SEQ # Path a donde colocar la secuencia filtrada (.fasta)
# SPLICE # Tener en cuenta el aplicing en el alineamiento (Bool) [true]
# EXONS # Extraer todo el rango de la secuencia alineada (Bool) [true]
# GENEXC # Excluir el cuerpo del gen [false]
# UTRS # Extraer las UTRs en la secuencia extraida (Bool) [false]
# UPBASES # Cuantas bases extraer upstream 5' [default 200]
# DOWNBASES # Cuantas bases extraer downstream 3' [default 200]

# Alingment of query to reference
if [ $ALN = false ] ; then
	echo "Query provided ($QUERY). Starting alignment to reference..."
	echo "Alignment command used:"
	echo "magicblast -splice $SPLICE -query $QUERY -subject $REF -out $OUTALN -no_unaligned"
        magicblast -splice $SPLICE -query $QUERY -subject $REF -out $OUTALN -no_unaligned && echo "Alignment successful"
	if [ $? -gt 0 ]; then echo "Problems with alignment, exiting..."; exit 1 ; fi
	echo "Alignment file is in $OUTALN"
else ## If an alingment is provided skip the process
	echo "Alignment provided ($ALN). Skipping alignment process."
        OUTALN=$ALN
fi

if [ $EXONS = "on" ] ; then ## In this case we only extract the exons of the gene
	echo "Extracting exons..."
	bamToBed -split -i $OUTALN | bedtools getfasta -s -fi $REF -bed - > $SEQ
	echo "Done. Extracted sequences are in $SEQ"
else
	if [ $UTRS = "on" ] ; then ## In this case we will extract the UTRs
		echo "Taking $UPBASES upstream and $DOWNBASES downstream bases."
		if [ $GENEXC = flank ]; then echo "Extracting only UTRs... "; else echo "Extracting the whole gene with UTRs..."; fi
		bamToBed -i $OUTALN | bedtools $GENEXC -s -l $UPBASES -r $DOWNBASES -i - -g <(bioawk -c fastx '{printf("%s\t%s\n", $name, length($seq))}' $REF) | bedtools getfasta -s -fi $REF -bed - > $SEQ
#		aux=`bioawk -c fastx '{printf("%s\t%s\n", $name, length($seq))}' $REF`
#		bamToBed -i $OUTALN | bedtools $GENEXC -s -l $UPBASES -r $DOWNBASES -i - -g $aux | bedtools getfasta -fi $REF -bed - > $SEQ
	else ## In this case we extract the whole gene alone
		echo "Extracting the whole gene..."
		bamToBed -i $OUTALN | bedtools getfasta -s -fi $REF -bed - > $SEQ
		echo "Done. Extracted sequences are in $SEQ"
	fi
fi
exit 0


