function [AssetDemand, Distribution] = LSDemand(r, Value, param, bound)
    AGridN = size(Value, 1);
    ZGridN = size(Value, 2);
    beta = param(1);
    [Value, UtilityMatrix, MarkovMatrix, AssetGrid] = LSHuggett(param, r, AGridN, ZGridN, Value, bound);
    [Policy, Index] = PolicyHuggett(beta, Value, AssetGrid,ZGridN,UtilityMatrix,MarkovMatrix);
    
    Distribution = InvariantDistribution(Index, MarkovMatrix);

    

    AssetDemand = sum(Policy.*Distribution, "all");
end
