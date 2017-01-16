scripts for getting FIMO pvalues/scores: 


Step 1)
# command for getting fasta sequences for all peaks
python get_fasta.pyImmGenATAC1219.peak.bed ImmGen.ATAC.master.peaks.fa

Step 2)
# command for running FIMO with CISBP database
fimo --text ~/motifDB/cisBP.meme ImmGen.ATAC.master.peaks.fa > smFIMO_cisBP_1219.txt

Step 3)
# command for creating a compact result file
cut -f1,2,6,7 smFIMO_cisBP_1219.txt > fimo_cisBP_scores.txt

Step 4)
# get min p-value per peak-motif combination
python get_min_pvalue_FIMO.py fimo_cisBP_scores.txt fimo_cisBP_scores_minPval.txt

Step 5)
# filter for Jason’s motif list
python filter_cisbp_jason.py fimo_cisBP_scores_minPval.txt fimo_cisBP_jmotifs.txt

Step 6)
# filter for associations with  p<1e-4
python threshold_fimo.py fimo_cisBP_jmotifs.txt 0.0001 fimo_cisBP_jmotifs_p4.txt
