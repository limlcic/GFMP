function [ outall ] = Run_conbine_steps( P,genesymbol )
%function [ out,Evaluation_value_all,Hub_module_P_all,Best_modules,Best_modules_gene ] = Run_conbine_steps( P,genesymbol )

outall.error = 0;
%%
%��һ�����ҳ�����ģ��
outall.member = P;
outall.symbols = genesymbol;

[ out ] = conbine_step1( P );%d�õ���out��ÿ�о���һ���໥������ģ��

if numel(out) == 0
    outall.error = 1;
    return;
end

%%
%��ʼ��outall
for i = 1:numel(out)
    outall.seed_module(i).pairs = out{i};
end
%%
%�ڶ������ϲ�
nm=size(P,2);
Raw_module=[1:1:nm];
for J=1:size(out,1)
    Hub_module=out{J,1};
    Rem_module=setdiff(Raw_module,Hub_module);
    Hub_module_P=P(:,Hub_module);
    Rem_module_P=P(:,Rem_module);
    [ Hub_module_P , history ] = conbine_step2 (Rem_module_P,Hub_module_P );
    Hub_module_P_all{J,1}=Hub_module_P;%�õ�ÿһ���໥����ģ�龭���ϲ����P����
    
    outall.seed_module(J).history = history;
    outall.seed_module(J).updated_member = Hub_module_P;
end

%%
%������������
AA=Hub_module_P_all;
for J=1:size(AA,1);
    A=AA{J,1};
    [ Evaluation_value ] = conbine_step3( A ) ;
    Evaluation_value_all{J,1}=Evaluation_value;%�õ����е��������
    
    outall.seed_module(J).eval_value = Evaluation_value;
end


%%
%���������ķ���ֵѰ������ģ��
Evaluation_value_all=cell2mat(Evaluation_value_all);
minus=min(Evaluation_value_all);
index=find(Evaluation_value_all==minus);
Best_modules=Hub_module_P_all(index,1);

Best_modules_gene=cell(size(Best_modules,1),1);
for JJ=1:size(Best_modules,1);
    %List=cell(size(Best_modules,1),1);
    List=cell(size(Best_modules{JJ,1},2),1);
    for II=1:size(Best_modules{JJ,1},2)
        aa=Best_modules{JJ,1}(:,II);
        bb=(find(aa>0));
       List{II} =genesymbol(bb);
    end
    Best_modules_gene{JJ}=List;%��ȡ�����ŷ�ģ�鷽���Ļ���ģ�顣
end


outall.best_modules = Best_modules;
% outall.best_modules_num = index;
outall.best_modules_num = out(index);
outall.best_modules_symbol = Best_modules_gene;

end

