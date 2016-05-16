function [Testp_sort, uniq_gene_expre_sort,uniq_gene_symbol_sort] = ttest_for_raw_data(gene_expre,gene_symbol,pn,p_val_thres)
%Run a simple t-test for the raw data
%
%Parameters:
%gene_expre:    m*n matrix of gene expression data 
%               m is the number of genes, n is the number of samples
%gene_symbol:   m*1 cell matrix of gene symbol
%pn             pn is the number of positive (negetive) samples
%p_val_thres:	The threshold of selected gene (default: 0)
%               If the value of threshold is larger than 1 (e.g. 
%               p_val_thres = 10), then the top 10 genes which have minimum
%               p values would be picked. If the value is between 0 and 1 
%               (e.g. p_val_thres = 0.05), then the genes whose p-values 
%               less than 0.05 would be picked. If the value is 0, then all
%               the gene would be picked.
%
%Output:
%uniq_gene_expre_sort:  Selected and sorted expression data
%uniq_gene_symbol_sort: Selected and sorted gene symbol




if nargin < 4
    p_val_thres = 0;
end



Testp={};
for i=1:size(gene_expre,1);
    [h,p]=ttest2(gene_expre(i,1:pn),gene_expre(i,pn+1:end));
    Testp=[Testp;p];
end
Testp=cell2mat(Testp);
[Testp_sort,bb]=sort(Testp);
uniq_gene_expre_sort=gene_expre(bb,:);
uniq_gene_symbol_sort=gene_symbol(bb,:);
if p_val_thres >= 1
    p_val_thres = round(p_val_thres);
    uniq_gene_expre_sort=uniq_gene_expre_sort(1:p_val_thres,:);
    uniq_gene_symbol_sort=uniq_gene_symbol_sort(1:p_val_thres,:);
elseif p_val_thres < 1 && p_val_thres > 0
    t = Testp_sort <= p_val_thres;
    uniq_gene_expre_sort = uniq_gene_expre_sort(t,:);
    uniq_gene_symbol_sort = uniq_gene_symbol_sort(t,:);
end



end

