function [ out ] = conbine_step1( member )
%合并的第一步，通过输入的矩阵确定好原始的模块分布
%member为m*n的矩阵，m为基因数，n为初始模块数，member中0元素表示该模块中没有这个基因
%输出为p*1的cell矩阵，每个cell为一个组合

%组合的判断标准是，每一组内两两不重叠，在这个基础上覆盖率尽量高
[m , n] = size(member);
out = cell(0);

%获取overlap
overlap = find_overlap(member);
overlap = overlap + overlap';

%初始化一个基础变量，用于结果输出
order = 1:n;

%开始统计模块
for i = 1:n
    sub_out = find_seed_m(overlap , i);
    out = [out;sub_out];
end

%判断重复
for i = 1:numel(out)
    out{i} = num2str(out{i});
end
out = unique(out);
for i = 1:numel(out)
    out{i} = str2num(out{i});
end

end

function out = find_seed_m(overlap , seed_line)
t = overlap(seed_line,:) == 0;
t(seed_line) = 0;
t = find(t == 1);
if numel(t) == 1
    out = {unique([seed_line,t])};
elseif numel(t) == 0
    out = cell(0);
else
    out = cell(0);
    %之前的算法太复杂，寻找种子节点我还是增强易读性好了
    %小范围内使用穷举法
    for i = numel(t): -1 : 1
        if i == 1
            group = t';
        else
            group = nchoosek(t,i);
        end
        for j = 1 : size(group,1)
            %还要确定是否存在包含关系
            c_handle = 0;
            for l = 1:numel(out)
                if numel(intersect( group(j,:) , out{l})) == numel(group(j,:))
                    c_handle = 1;
                    break;
                end
            end
            if c_handle
                continue;
            end
            if i == 1
                out = [out;{group(j,:)}];
                continue;
            end
            isok = 1;
            positions = nchoosek(group(j,:),2);
            for k = 1 : size(positions,1)
                if overlap( positions(k,1) , positions(k,2) ) ~= 0
                    isok = 0;
                end
            end
            if isok            
                out = [out;{group(j,:)}];
            end
        end
    end
    
    
    for i = 1:numel(out)
        out{i} = unique([out{i},seed_line]);
    end
end

end

function overlap = find_overlap(member)
[m , n] = size(member);
overlap = zeros(n , n);
for i = 1 : n
    for j = i + 1 : n
        t = member(: , i) .* member(: , j);
        s = member(: , i) + member(: , j);
        overlap(i , j) = sum(t ~= 0) / sum(s ~= 0);
%         s1 = sum(t ~= 0) / sum(member(: , i) ~= 0);
%         s2 = sum(t ~= 0) / sum(member(: , j) ~= 0);
%         overlap(i , j) = max(s1 , s2);
        
    end
end

end

function out = if_repeat(pairs , t)
out = 0;
if numel(pairs) < 1
    out = 0;
else
    for i = 1 : numel(pairs)
        s = pairs{i};
        out = out + sum(s - t);
    end
end

end


