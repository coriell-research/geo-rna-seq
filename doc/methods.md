## Methods

The Gene Expression Omnibus (GEO) [1]. was queried for BioProjects with terms 
related to 'cytotoxic drugs', 'chemotherapy', 'cancer', and 'rna-seq'. 
BioProjects with associated raw fastq files and metadata were downloaded from 
The Sequence Read Archive (SRA) [2]. For each BioProject, raw fastq files were 
uniformly processed using the REdiscoverTE pipeline [3]. Briefly, the REdiscoverTE
pipeline allows for the simultaneous quantification of gene and transposable 
element expression by aligning reads with Salmon [4] to a custom transcriptome 
consisting of all distinct RNA transcripts in GENCODE v26, RepeatMasker 
elements (n = 5,099,056) from the standard chromosomes, excluding all polyA 
repetitive elements, and distinct sequences representing GENCODE RE-containing 
introns (n = 185,403) and excluding any regions overlapping with exons on either 
strand. For each BioProject, gene and repetitive element counts were then normalized
using TMM normalization in edgeR [5] and differential expression testing was 
performed for each drug relative to it's control. A meta-analysis of all BioProjects
was then performed using the differential expression testing results for each study.  


## References

1.  Edgar R, Domrachev M, Lash AE. Gene Expression Omnibus: NCBI gene expression and hybridization array data repository Nucleic Acids Res. 2002 Jan 1;30(1):207-10
2. Leinonen R, Sugawara H, Shumway M; International Nucleotide Sequence Database Collaboration. The sequence read archive. Nucleic Acids Res. 2011;39(Database issue):D19-D21. doi:10.1093/nar/gkq1019
3. Kong, Y., Rose, C.M., Cass, A.A. et al. Transposable element expression in tumors is associated with immune infiltration and increased antigenicity. Nat Commun 10, 5228 (2019). https://doi.org/10.1038/s41467-019-13035-2
4. Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware quantification of transcript expression. Nat Methods. 2017;14(4):417-419. doi:10.1038/nmeth.4197
5. Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics, 26(1), 139-140. doi: 10.1093/bioinformatics/btp616. 