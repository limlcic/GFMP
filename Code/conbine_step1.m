function [ out ] = conbine_step1( member )
%�ϲ��ĵ�һ����ͨ������ľ���ȷ����ԭʼ��ģ��ֲ�
%memberΪm*n�ľ���mΪ��������nΪ��ʼģ������member��0Ԫ�ر�ʾ��ģ����û���������
%���Ϊp*1��cell����ÿ��cellΪһ�����

%��ϵ��жϱ�׼�ǣ�ÿһ�����������ص�������������ϸ����ʾ�����
[m , n] = size(member);
out = cell(0);

%��ȡoverlap
overlap = find_overlap(member);
overlap = overlap + overlap';

%��ʼ��һ���������������ڽ�����
order = 1:n;

%��ʼͳ��ģ��
for i = 1:n
    sub_out = find_seed_m(overlap , i);
    out = [out;sub_out];
end

%�ж��ظ�
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
    %֮ǰ���㷨̫���ӣ�Ѱ�����ӽڵ��һ�����ǿ�׶��Ժ���
    %С��Χ��ʹ����ٷ�
    for i = numel(t): -1 : 1
        if i == 1
            group = t';
        else
            group = nchoosek(t,i);
        end
        for j = 1 : size(group,1)
            %��Ҫȷ���Ƿ���ڰ�����ϵ
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


