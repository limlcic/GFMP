function [ top20_P,top20_hubs ] = get_top_P( genesymbol,C,index1, topNum )
%找出一圈度，排名前20的基因所在的模块
%  genesymbopl是对应的原来条全出的差异性表达基因，C是上一步得到的归一化矩阵，all_P是得到的所有基因所在的P矩阵，all_module对应着每一个基因所在的模块。index1是度排序后对应的位置
all_P={};
for I=1:size(genesymbol,1);
M=zeros(size(genesymbol,1),1);%M为权重矩阵。
M(I,1)=1;
P1=C*M;
nonzero=find(P1~=0);
all_module{I,1}=genesymbol(nonzero);
all_P=[all_P;P1];
end
top20_modules=all_module(index1(1:topNum,1));
top20_hubs=genesymbol(index1(1:topNum,1));
top20_P=all_P(index1(1:topNum,1));
end

