function [ Hub_module_P , history] = conbine_step2 (Rem_module_P,Hub_module_P )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    for j=1:size(Rem_module_P,2);
       for i=1:size(Hub_module_P,2);
            %nr=numel(find(Rem_module_P(:,j)~=0));
            cc=Rem_module_P(:,j).*Hub_module_P(:,i);
            nc=numel(find(cc~=0));%交集的基因个数
            uu=Rem_module_P(:,j)+Hub_module_P(:,i);
            nu=numel(find(uu~=0));%并集的基因个数
            nn=nc/nu;
            overlap(i,j)=nn;%overlap一共是i行*j列，每一列代表其他模块与Hub_module的重叠率
       end
       %idx0=find(overlap(:,j)==max(overlap(:,j)));%判断每一列数据中最大的值的位置就是要和Hub_module的位置
      % Hub_module_P(:,idx0)=Rem_module_P(:,j)+Hub_module_P(:,idx0);%合并模块，然后下一个模块
        idx=find(overlap(:,j)==max(overlap(:,j)));%判断每一列数据中最大的值的位置就是要和Hub_module的位置
        idx0=idx(1,1);
       Hub_module_P(:,idx0)=Rem_module_P(:,j)+Hub_module_P(:,idx0);%合并模块，然后下一个模块
       history(j).overlap = overlap(:,end);
       history(j).selected = idx0;
    end    
end



