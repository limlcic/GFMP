function [ Evaluation_value ] = conbine_step3( A )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
for i=1:size(A,2);
    for j=1:size(A,2);
        commom_gene=A(:,i).*A(:,j);%���о������
        n0=numel(find(commom_gene~=0));%������ֵ��������ģ��֮�乲�еĻ���
        n1=numel(find(A(:,i)~=0));%����һ��ģ��Ļ������
        n2=numel(find(A(:,j)~=0));%����һ������ģ��ĸ���
        n=min(n1,n2);%����ģ���л������С��ģ��������
    overlap{i,j}=n0/n;%�õ��ϴ�ֵ��overlap���Խ�����Ϊ1 
    end
end
%Evaluation_value=min(min(cell2mat(overlap)));%�ҳ���Խ�С��overlapֵ��Ϊ������׼��
N=size(A,2);
T=N*(N-1)/2;
M=triu(cell2mat(overlap),1);%ȡ�϶Խ��ߵ�Ԫ��(�������Խ���)�����ǰ��ص��ʶ�ȡ����
Evaluation_value= sum(sum(M))/T;
end

