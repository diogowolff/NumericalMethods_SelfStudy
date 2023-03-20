function [r, Value, Policy] = HuggettEquilibrium(param, AGridN, ZGridN)
    beta = param(1);

    RLowerBound = .0;
    RUpperBound = 1/beta-1;
    
    error = 1;
    tol = 10e-6;
    iter = 0;

    while error >= tol & iter < 500
        r = (RLowerBound + RUpperBound)/2;
        
        if iter>0
            [Value, UtilityMatrix, MarkovMatrix, AssetGrid] = ExperimentalHuggett(param, r, AGridN, ZGridN, Value);
        else
            [Value, UtilityMatrix, MarkovMatrix, AssetGrid] = ExperimentalHuggett(param, r, AGridN, ZGridN);
        end

        [Policy, Index] = PolicyHuggett(beta, Value, AssetGrid,ZGridN,UtilityMatrix,MarkovMatrix);
    
        Distribution = InvariantDistribution(Index, MarkovMatrix);
    
        AssetDemand = sum(Policy.*Distribution, "all");
        
        if AssetDemand>0
            RUpperBound = r;
        else
            RLowerBound = r;
        end

        error = abs(r);
        display([AssetDemand, RLowerBound, RUpperBound]);
        iter = iter+1;
    end

end

