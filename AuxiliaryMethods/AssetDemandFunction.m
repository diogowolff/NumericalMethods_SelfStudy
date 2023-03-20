function [AssetDemand, Distribution] = AssetDemandFunction(r, Value, Distribution)
    param = [.96, 1.0001, .9, .01];
    AGridN = size(Value, 1);
    ZGridN = size(Value, 2);
    beta = param(1);
    [Value, UtilityMatrix, MarkovMatrix, AssetGrid] = ExperimentalHuggett(param, r, AGridN, ZGridN, Value);
    [Policy, Index] = PolicyHuggett(beta, Value, AssetGrid,ZGridN,UtilityMatrix,MarkovMatrix);
    
        Distribution = InvariantDistribution(Index, MarkovMatrix);

    

    AssetDemand = sum(Policy.*Distribution, "all");
end

mu = (eye(size(TransitionMatrix, 1)) - TransitionMatrix + ones(size(TransitionMatrix)))\ones(size(TransitionMatrix), 1);