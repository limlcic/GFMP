function [ outall ] = divide_network( theta,genesymbol,topNum  )
%Divede the constructed network into modules
%Parameters:
%theta: m*m netowrk matrix (adjacent matrix). m is the number of genes.
%genesymbol:    m*1 cell matrix of gene symbol
%topNum:        The number of picked hub genes (default: 20)
%               If topNum = 20, then the first 20 gene which have maximum
%               degree would be picked for building seed modules.
% 
%Output:        
%out:               A structure of the selected modules
%   out.error:      If the value is 1, then means the data can't be used
%                   for module partition.
%   out.member:     m*k matrix. m is the number of gene, k is the number of
%                   seed modules. The genes in the modules are represented 
%                   as nonzero element in the matrix For example,
%                   out.member(45,3) nonzero means the 45th gene is in the
%                   third module.
%   out.symbols:    The symbols of the genes
%       out.seed_module:    The seed modules and the related parameters
%       out.seed_module.pairs:  The number is related with the row in
%                               out.member. The modules are the begin of combine.
%       out.seed_module.history:    The detailed step of combine.
%       out.seed_module.updated_member: The nonzero value means belonging. The
%                                       updated_member is similar with the
%                                       out.member.
%       out.seed_module.eval_value: The evaluation result of the combined
%                                   modules.
%   out.best_modules: The modules which have minimum eval_value.
%   out.best_modules_num:   The related pairs of the best modules.
%   out.best_modules_symbol:    The related symbols of the modules.

if nargin < 3
    topNum = 20;
    if nargin < 2
        genesymbol = num2cell((1:size(theta,1))');
        for i = 1:numel(genesymbol)
            genesymbol{i} = num2str(genesymbol{i});
        end
    else
        genesymbol = importdata(genesymbol);
    end
end

[ C ] = get_C_vetor(theta);
[ new_degree ] = get_new_degree( theta );
[rank21,index1]=sort(new_degree,'descend');
[ top20_P ] = get_top_P( genesymbol,C,index1,topNum );
P=top20_P;
P=P';
P=cell2mat(P);
[ outall] = Run_conbine_steps( P,genesymbol );

end

