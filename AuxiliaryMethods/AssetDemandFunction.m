function [AssetDemand] = AssetDemandFunction(r, Value)
    param = [.96, 1.0001, .9, .01];
    AGridN = 100;
    ZGridN = 9;
    beta = param(1);
    [Value, UtilityMatrix, MarkovMatrix, AssetGrid] = HowardHuggett(param, r, AGridN, ZGridN, Value);
    [Policy, Index] = PolicyHuggett(beta, Value, AssetGrid,ZGridN,UtilityMatrix,MarkovMatrix);
    
    Distribution = InvariantDistribution(Index, MarkovMatrix);

    AssetDemand = sum(Policy.*Distribution, "all");
end

