function [ top20_P,top20_hubs ] = get_top_P( genesymbol,C,index1, topNum )
%�ҳ�һȦ�ȣ�����ǰ20�Ļ������ڵ�ģ��
%  genesymbopl�Ƕ�Ӧ��ԭ����ȫ���Ĳ����Ա�����C����һ���õ��Ĺ�һ������all_P�ǵõ������л������ڵ�P����all_module��Ӧ��ÿһ���������ڵ�ģ�顣index1�Ƕ�������Ӧ��λ��
all_P={};
for I=1:size(genesymbol,1);
M=zeros(size(genesymbol,1),1);%MΪȨ�ؾ���
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

