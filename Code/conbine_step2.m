function [ Hub_module_P , history] = conbine_step2 (Rem_module_P,Hub_module_P )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    for j=1:size(Rem_module_P,2);
       for i=1:size(Hub_module_P,2);
            %nr=numel(find(Rem_module_P(:,j)~=0));
            cc=Rem_module_P(:,j).*Hub_module_P(:,i);
            nc=numel(find(cc~=0));%�����Ļ������
            uu=Rem_module_P(:,j)+Hub_module_P(:,i);
            nu=numel(find(uu~=0));%�����Ļ������
            nn=nc/nu;
            overlap(i,j)=nn;%overlapһ����i��*j�У�ÿһ�д�������ģ����Hub_module���ص���
       end
       %idx0=find(overlap(:,j)==max(overlap(:,j)));%�ж�ÿһ������������ֵ��λ�þ���Ҫ��Hub_module��λ��
      % Hub_module_P(:,idx0)=Rem_module_P(:,j)+Hub_module_P(:,idx0);%�ϲ�ģ�飬Ȼ����һ��ģ��
        idx=find(overlap(:,j)==max(overlap(:,j)));%�ж�ÿһ������������ֵ��λ�þ���Ҫ��Hub_module��λ��
        idx0=idx(1,1);
       Hub_module_P(:,idx0)=Rem_module_P(:,j)+Hub_module_P(:,idx0);%�ϲ�ģ�飬Ȼ����һ��ģ��
       history(j).overlap = overlap(:,end);
       history(j).selected = idx0;
    end    
end



