function [ new_degree ] = get_new_degree( T )
%T是graphicallasso计算出来的相关系数矩阵
%  这里得到的度是一阶的，并且去掉了自相关
theta=T;
new_theta=theta-diag(diag(theta));
new_network=sparse(new_theta);
[rowId, colId]=find(new_network~=0);
new_data=[rowId,colId,full(new_network(new_network~=0))];
new_network1=full(new_network);
new_network1(new_network1~=0)=1; 
new_degree=sum(new_network1,2);
end

