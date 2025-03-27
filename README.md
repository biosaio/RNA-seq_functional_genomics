**RNA-seq project**

In this project we were tasked with choosing a functional genomics RNA-seq dataset, and to process it using deseq2 pipeline. We were further required to derive useful insights from the resulting analysis and present
our findings to the professor and fellow students. 

Here I include all R files and the final knit with the project details, results and their discussion.

- setup.R : in this file we make sure that all packages that will be used in the project are properly installed and loaded
- GEOquery.R : Downloads the dataset from Gene Expression Omnibus
- DESeq2DataSet.R : Creates the DESeq2 dataset starting from the downloaded RNA seq data
- DESeqDiffExprAnal.R : Here we run the DEseq2 analysis and obtain the first visualizations
- DESeq2visualization.R : use log transformation to reduce variance in low conunt genes. Plot results
- Annotation.R : annotate genes with ensembl IDs. Volcano plot to display differentially expressed genes
- GoGeneSetEnriAnal.R : GO enrichment analysis. Visualize which GO are differentially expressed
- KeggGeneSetEnriAnal.R : Investigate which KEGG patways are impacted by differentially expressed genes

To have a complete overview of the project, I suggest consulting FGProjectRMarkdown, which will go through all the relevant code and give thorough explanations and interpretation for the results
