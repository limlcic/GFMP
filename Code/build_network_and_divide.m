function [ out ] = build_network_and_divide( genesymbol,data,LambdaStart , LambdaEnd, StepSize, topNum)

%Build a netowrk by using graphic lasso and then divide it into modules
%
%Parameters:
%gene_symbol:   m*1 cell matrix of gene symbol
%data:          m*n matrix of gene expression data 
%               m is the number of genes, n is the number of samples
%
%LambdaStart LambdaEnd StepSize:
%               The Lambda is a parameter for graphic lasso, the three
%               parameters is to find the best value of lambda. If
%               LambdaStart = 0.1, LambdaEnd = 1 and StepSize = 0.2, then
%               5 networks with different Lambda (0.1 0.3 ... 0.9) would
%               be built.
%
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
%

if nargin < 6
    topNum = 20;
        if nargin < 5
            StepSize = 0.1;
            if nargin < 4
                LambdaEnd = 1;
                if nargin < 3
                    LambdaStart = 0.1;
                end
            end
        end
end




% data = data';
% genesymbol = genesymbol';

for lambda = LambdaStart : StepSize : LambdaEnd
    [w, theta, iter, avgTol, hasError] = GraphicalLasso(data', lambda);
    [ C ] = get_C_vetor(theta);
    [ new_degree ] = get_new_degree( theta );
    [rank21,index1]=sort(new_degree,'descend');
    [ top20_P ] = get_top_P( genesymbol,C,index1, topNum );
    P=top20_P;
    P=P';
    P=cell2mat(P);
    [ outall] = Run_conbine_steps( P,genesymbol );
    if outall.error ~= 1
        out = outall;
        out.lambda = lambda;
        out.theta = theta;
        return;
    end
end
out = outall;


end

