function [Policy, Index] = PolicyHuggett(beta, Value, AssetGrid,ZGridN,UtilityMatrix,MarkovMatrix)

    AGridN = size(AssetGrid, 2);
    Index = zeros(AGridN, ZGridN);
    Policy = zeros(AGridN, ZGridN);
    
    for i=1:AGridN
        for j=1:ZGridN
            PossibleValues = UtilityMatrix(:, i, j) + ...
                    beta.*Value*MarkovMatrix(j,:)';
            [~, Index(i, j)] = max(PossibleValues,[],'all',"linear");
            Policy(i,j) = AssetGrid(Index(i, j));
        end
    end
end

