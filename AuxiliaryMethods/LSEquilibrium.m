function [r, Value, Policy] = LSEquilibrium(param, AGridN, ZGridN, bound)
    beta = param(1);

    RLowerBound = .0;
    RUpperBound = 1/beta-1;
    
    error = 1;
    tol = 10e-6;
    iter = 0;

    while error >= tol & iter < 50
        r = (RLowerBound + RUpperBound)/2;
        
        if iter>0
            [Value, UtilityMatrix, MarkovMatrix, AssetGrid] = LSHuggett(param, r, AGridN, ZGridN, Value, bound);
        else
            [Value, UtilityMatrix, MarkovMatrix, AssetGrid] = LSHuggett(param, r, AGridN, ZGridN, bound);
        end

        [Policy, Index] = PolicyHuggett(beta, Value, AssetGrid,ZGridN,UtilityMatrix,MarkovMatrix);
    
        Distribution = InvariantDistribution(Index, MarkovMatrix);
    
        AssetDemand = sum(Policy.*Distribution, "all");
        
        if AssetDemand>0
            RUpperBound = r;
        else
            RLowerBound = r;
        end

        error = abs(RUpperBound - RLowerBound);
        display([AssetDemand, RLowerBound, RUpperBound]);
        iter = iter+1;
    end

end
