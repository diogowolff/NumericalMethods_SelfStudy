function [InitialGuess] = InvariantDistribution(PolicyIndex, MarkovMatrix)
    DistDims = size(PolicyIndex);
    InitialGuess = ones(DistDims)./numel(PolicyIndex);

    
    NextGuess = InitialGuess;

    error = 1;
    tol = 10e-10;
    iter = 0;

    while error >= tol & iter < 1000
        for i=1:size(PolicyIndex,1)
            for j=1:size(PolicyIndex,2)
                NextGuess(i,j) = sum(InitialGuess.*(PolicyIndex == i))*MarkovMatrix(:,j);
            end
        end
        
        error = max(abs(InitialGuess - NextGuess),[],'all',"linear");
        %disp(error);
        InitialGuess = NextGuess;
        iter = iter+1;
    end
end

