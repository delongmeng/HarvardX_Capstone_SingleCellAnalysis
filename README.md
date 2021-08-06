# Single Cell RNA and Protein Profiling Using Seurat

Can also be visited here:
http://rpubs.com/delongmeng/563428

# Introduction

Single cell RNA-sequencing (scRNA-seq) is a new technology that is rapidly growing recently in the biomedical science field. Prior to scRNA-seq, people use so-called bulk RNA-seq to perform transcriptome analysis, namely levels of gene expression of the whole genome in different groups of biomedical samples, such as tumor samples and normal control tissues, or cells of different treatments. By comparing gene expression in different samples, we usually get a list of differentially expressed genes, associated with disease or cell status. However, this bulk RNA-seq technique ignores the heterogeneity within a group and only capture the average gene expression of different cell in a group.

Through scRNA-seq, we can obtain gene expression levels of the whole genome in each individual cell. This is achived by labeling each cell with a unique barcode of DNA sequence. When we amplify the complementary DNA (cDNA) which is synthesized from RNAs, the unique barcode will also be amplified and eventually provide the cell identity information after the amplified DNAs are sequenced. The task of scRNA-seq data analysis is to perform an unsupervised machine learning process, and explore and identify the patterns of variation in cell states. A common strategy is to reduce the dimentionality and cluster the cells into discrete groups in a 2-D space, and identify the corresponding gene expression markers that distinguish these separated groups. We believe that cell populations with similar molecular features might also have similar functions and behaviors under pathophysiological conditions. 

For this Capstone project, I will use a public dataset of scRNA-seq study of human cord blood mononuclear cells (CBMCs). Interestingly, the authors of this study creatively developed a new approach to simultaneously analyze RNA expression levels and cell surface protein levels (i.e., epitopes), which they called "cellular indexing of transcriptomes and epitopes by sequencing" (CITE-seq) [1]. Cell surface protein markers are critical indicators of immunophenotypes of blood cells, and are usually measured by flow cytometry, in which fluorescence labeled antibodies will be used to recognize these protein markers and quantified based on the fluorescence. The brief concept of this CITE-seq is that they conjugate the antibodies to oligonucleotides (oligos) that can be captured in the process of scRNA-seq library preparations. This oligo contains a unique barcode that can provide the information of the antibody in later analysis, and can be released from the antibody during cell lysis, generating so-called antibody-derived tags (ADTs) that can be labeled by the same cellular barcode with the RNAs and amplified together with the RNAs. This new approach makes it possible to analyze different levels of molecular features (RNAs and proteins here) at the same time, and achieve more detailed characterization of cellular phenotypes.

Here I will mainly use the Seurat R package to perform the single cell analysis. Seurat is a popular and powerful toolkit that provides solutions for each stage of single cell analysis. Developed by the *Satija lab*, now the most recent version is *Seurat Version 3* [2]. Very helpful tutorials can also be found in their webpage: https://satijalab.org/seurat/.

# References

[1] Stoeckius, M., Hafemeister, C., Stephenson, W. et al. Simultaneous epitope and transcriptome measurement in single cells. Nat Methods 14, 865–868 (2017) doi:10.1038/nmeth.4380

[2] Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E., Mauck III, W. M. et al. Comprehensive Integration of Single-Cell Data. Cell 177, 1888-1902 (2019) doi:10.1016/j.cell.2019.05.031. 
