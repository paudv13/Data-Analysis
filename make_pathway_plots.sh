#!/bin/bash


#####USAGE#####
# ./make_pathway_plots.sh pathways_list_file.txt inputfile.tsv variable OUTPUT_DIRECTORY

pathways="/mnt/synology/GERARD/DATA/PhD/HUMANN_TABLES/da_pwy_cd_hc.list"
infile="/mnt/synology/GERARD/DATA/PhD/HUMANN_TABLES/rna_pwys_humann.barplot.tsv"
variable="DIS"
outdir="/home/gsergom/PATHWAY_PLOTS"

pathways=$1
infile=$2
variable=$3
outdir=$4

mkdir -v $outdir

while read p
do
  humann_barplot --input ${infile} --focal-feature ${p} --focal-metadata ${variable} --last-metadata ${variable} \
  --output ${outdir}"/barplot_"${p}".png" --sort sum metadata
  MYVARIABLE="$(humann_barplot --input ${infile} --focal-feature ${p} --focal-metadata ${variable} --last-metadata ${variable} \
  --output ${outdir}"/barplot_"${p}".png" --sort sum metadata 2>&1 > /dev/null)"
#  humann_barplot --input ${infile} --focal-feature ${p} --focal-metadatum ${variable} --last-metadatum ${variable} --output ${outdir}"/barplot_"${p}".png"
#  humann2_barplot --input infile --focal-feature ${p} --focal-metadatum variable --last-metadatum variabl --output "barplot_"${p}".png"
  if [[ $MYVARIABLE != "" ]]; then
    echo "Could not plot ${p}, this might be because there is no species information about the pathway..."
  else
    echo "Plotting pathway ${p}..."
  fi
done <$pathways
