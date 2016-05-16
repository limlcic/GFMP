function [out] = Get_result( genesymbol,data,LambdaStart , LambdaEnd, StepSize)

if nargin < 5
    StepSize = 0.1;
    if nargin < 4
        LambdaEnd = 1;
        if nargin < 3
            LambdaStart = 0.1;
        end
    end
end


for lambda = LambdaStart : StepSize : LambdaEnd
    [w, theta, iter, avgTol, hasError] = GraphicalLasso(data, lambda);
    [ C ] = get_C_vetor(theta);
    [ new_degree ] = get_new_degree( theta );
    [rank21,index1]=sort(new_degree,'descend');
    [ top20_P ] = get_top_P( genesymbol,C,index1 );
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

end

