function [ C ] = get_C_vetor(T )
A1=T;
B1=abs(A1);
[m n]=size(B1);
CC={};
for i=1:m;
  CC=[CC,B1(:,i)/sum(B1(:,i))];
end
C=cell2mat(CC);
end

