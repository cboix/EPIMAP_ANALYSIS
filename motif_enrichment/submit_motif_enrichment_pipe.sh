#!/bin/bash
# ---------------------------------
# Submit motif enrichment pipeline:
# ---------------------------------
source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh 18

# Variables and files:
MOTIFDIR=$DBDIR/motifs
ALL_PFM=$MOTIFDIR/collated_pfm.txt
MNUM=$( grep ">" $ALL_PFM | wc -l )

conda activate base

# Run motif matches:
qsub -cwd -t 1-$MNUM -l h_vmem=20G -l h_rt=1:00:00 -N motif_matches -j y -b y -V -r y -o $DBDIR/out/ "$BINDIR/motif_enrichment/run_motif_matches.sh -d $ALL_PFM -k 8"

# Run motif matches (hg38):
qsub -cwd -t 1-$MNUM -l h_vmem=20G -l h_rt=4:00:00 -N motif_matches -j y -b y -V -r y -o $DBDIR/out/ "$BINDIR/motif_enrichment/run_motif_matches.sh -d $ALL_PFM -k 8 -g hg38"

# Command to intersect with the mapping for each DHS:
# qsub -cwd -t 1-$MNUM -l h_vmem=30G -l h_rt=0:18:00 -hold_jid motif_matches -N reduce_matches -j y -b y -V -r y -o $DBDIR/out/ "$BINDIR/motif_enrichment/reduce_motif_matches.sh -d $ALL_PFM -k 8"
# qsub -cwd -t 1-$MNUM -l h_vmem=12G -l h_rt=0:20:00 -hold_jid motif_matches -N reduce_matches -j y -b y -V -r y -o $DBDIR/out/ "$BINDIR/motif_enrichment/reduce_motif_matches.sh -d $ALL_PFM -k 8"
qsub -cwd -t 1-$MNUM -l h_vmem=30G -l h_rt=0:45:00 -hold_jid motif_matches -N reduce_matches -j y -b y -V -r y -o $DBDIR/out/ "$BINDIR/motif_enrichment/reduce_motif_matches.sh -d $ALL_PFM -k 8"

# TODO: Collate all hits: take only ones with some threshold of either raw or ctrl pval and fold enrichment.

cd $MOTIFDIR
CLSMATCH=$MOTIFDIR/collated.bkgdhs.cls.enrich.tsv
rm $CLSMATCH $CLSMATCH.gz
for file in `ls matches/*.bkgdhs.cls.enrich.tsv.gz`; do
    echo $file
    if [[ ! -s $CLSMATCH ]]; then
        zcat $file >> $CLSMATCH
    else
        zcat $file | awk 'NR >1' >> $CLSMATCH
    fi
done

gzip $CLSMATCH


cd $MOTIFDIR
CLSMATCH=$MOTIFDIR/collated.bkgdhs.epi.enrich.tsv
rm $CLSMATCH $CLSMATCH.gz
for file in `ls matches/*.bkgdhs.epi.enrich.tsv.gz`; do
    echo $file
    if [[ ! -s $CLSMATCH ]]; then
        zcat $file >> $CLSMATCH
    else
        zcat $file | awk 'NR >1' >> $CLSMATCH
    fi
done

gzip $CLSMATCH


# Run missed motif matches:
if [[ "0" == "1" ]]; then
    while read ITER; do
        echo $ITER
        qsub -cwd -t $ITER -l h_vmem=40G -l h_rt=8:00:00 -N motif_matches_${ITER} -j y -b y -V -r y -o $DBDIR/out/ "$BINDIR/motif_enrichment/run_motif_matches.sh -d $ALL_PFM -k 8 -g hg38"
    done < rerun_matches.tsv

    # Ensure all run ok.
    while read ITER; do
        qsub -cwd -t $ITER -l h_vmem=40G -l h_rt=2:00:00 -hold_jid motif_matches_${ITER} -N reduce_matches_${ITER} -j y -b y -V -r y -o $DBDIR/out/ "$BINDIR/motif_enrichment/reduce_motif_matches.sh -d $ALL_PFM -k 8"
    done < rerun_matches.tsv
fi

