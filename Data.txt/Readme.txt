AML_gene_expression.txt is the AML raw data we used in this study, the row is the gene expression, the column is the sample
GBM_gene_expression.txt is the GBM(glioblastoma )raw data we used in this study, the row is the gene expression, the column is the sample
NB_gene_expression.txt is the NB(neuroblastoma )raw data we used in this study, the row is the gene expression, the column is the sample
NB_valid_gene_expression.txt is the independent testing set NB(neuroblastoma )raw data we used in this study, the row is the gene expression, the column  is the sample
GBM_seq_gene_expression.txt is the GBM_seq(the mRNA dataset of the glioblastoma)raw data we used in this study, the row is the gene expression, the column is the sample

AML_gene_symbol.txt is the AML gene that used in the study, row is the gene, which are corresponding to the row of gene expression 
GBM_gene_symbol.txt is the GBM gene that used in the study, row is the gene, which are corresponding to the row of gene expression 
NB_gene_symbol.txt is the NB gene that used in the study, row is the gene, which are corresponding to the row of gene expression 
NB_valid_gene_symbol.txt is the independent testing set NB gene that used in the study, row is the gene, which are corresponding to the row of gene expression 
GBM_seq_gene_symbol.txt is the GBM_seq gene that used in the study, row is the gene, which are corresponding to the row of gene expression 

AML_DEG_expression.txt is the AML DEG (differently expression genes) with ttest p-value <=0.05 that we selected from the gene expression
GBM_DEG_expression.txt is the GBM DEG (differently expression genes) with ttest p-value <=0.05 that we selected from the gene expression
NB_DEG_expression.txt is the NB DEG (differently expression genes) with ttest p-value <=0.05 that we selected from the gene expression
GBM_seq_DEG_expression.txt is the GBM_seq DEG (differently expression genes) with ttest p-value <=0.05 that we selected from the gene expression

AML_DEG_symbol.txt is the AML DEG (differently expression genes) that are corresponding to the DEG expression 
GBM_DEG_symbol.txt is the GBM DEG (differently expression genes) that are corresponding to the DEG expression 
NB_DEG_symbol.txt is the NB DEG (differently expression genes) that are corresponding to the DEG expression 
GBM_seq _DEG_symbol.txt is the GBM_seq DEG (differently expression genes) that are corresponding to the DEG expression 

AML_gene_network.txt is the AML DEG (differently expression genes) gene networks we constructed by using graphical lasso,the row are corresponding to the AML_DEG_symbol and the column are the same with row
GBM_gene_network.txt is the GBM DEG (differently expression genes) gene networks we constructed by using graphical lasso,the row are corresponding to the GBM_DEG_symbol and the column are the same with row
NB_gene_network.txt is the NB DEG (differently expression genes) gene networks we constructed by using graphical lasso,the row are corresponding to the NB_DEG_symbol and the column are the same with row
GBM_seq _gene_network.txt is the GBM_seq DEG (differently expression genes) gene networks we constructed by using graphical lasso, the row are corresponding to the NB_DEG_symbol and the column are the same with row

DEGs are selected from gene expression profiles.
DEGs are directly used to construct gene networks
Gene networks can be directly used to partition functional modules by using Functional Modules Partition algorithm