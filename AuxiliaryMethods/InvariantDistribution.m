function [DistDims,InitialGuess,NextGuess] = InvariantDistribution(PolicyIndex, MarkovMatrix)
    DistDims = size(PolicyIndex);
    InitialGuess = ones(DistDims)./numel(DistDims);
    NextGuess = InitialGuess;

    error = 1;
    tol = 10e-5;
    iter = 0;

    while error >= tol & iter <10
        for i=1:size(PolicyIndex,1)
            for j=1:size(PolicyIndex,2)
                NextGuess(i,j) = sum(InitialGuess.*(PolicyIndex == i))*MarkovMatrix(:,j);
            end
        end
        
        error = max(abs(InitialGuess - NextGuess),[],'all',"linear");
        disp(error);
        InitialGuess = NextGuess./sum(NextGuess, "all");
        iter = iter+1;
    end
end

