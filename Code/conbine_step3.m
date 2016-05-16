function [ Evaluation_value ] = conbine_step3( A )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
for i=1:size(A,2);
    for j=1:size(A,2);
        commom_gene=A(:,i).*A(:,j);%两列矩阵相乘
        n0=numel(find(commom_gene~=0));%非零数值代表两两模块之间共有的基因
        n1=numel(find(A(:,i)~=0));%其中一个模块的基因个数
        n2=numel(find(A(:,j)~=0));%另外一个基因模块的个数
        n=min(n1,n2);%两个模块中基因书较小的模块基因个数
    overlap{i,j}=n0/n;%得到较大值的overlap，对角线上为1 
    end
end
%Evaluation_value=min(min(cell2mat(overlap)));%找出相对较小的overlap值作为评估标准。
N=size(A,2);
T=N*(N-1)/2;
M=triu(cell2mat(overlap),1);%取上对角线的元素(不包括对角线)，就是吧重叠率都取出来
Evaluation_value= sum(sum(M))/T;
end

