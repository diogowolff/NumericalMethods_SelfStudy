function [Value, UtilityMatrix, MarkovMatrix, AssetGrid] = ExperimentalHuggett(param, r, AGridN, ZGridN, StartingGuess)
    
    beta = param(1);
    gamma = param(2);
    rho = param(3);
    sigma = param(4);
    
    

    [ShockGrid, MarkovMatrix] = TauchenDiscretizer(ZGridN, 3, 0, rho, sigma);
    
    phi = (exp(ShockGrid(1))/(1-beta)); 

    AssetGrid = linspace(-phi, phi, AGridN);
    [AGrid, ANewGrid, TFPGrid] = meshgrid(AssetGrid, AssetGrid, exp(ShockGrid));

    UtilityMatrix = ((TFPGrid + (1+r).*AGrid - ANewGrid).^(1-gamma)-1)./(1-gamma);
    UtilityMatrix(TFPGrid + (1+r).*AGrid - ANewGrid <=0 ) = -inf;

    if ~exist('StartingGuess','var')
        Value = zeros(AGridN, ZGridN);
    else
        Value = StartingGuess;
    end

    

    NewValue = Value;
    error = 100;
    tol = 10e-6;
    iter = 0;
    Index = Value;


    
    while error > tol
        if iter<50
            NewValue = reshape(max(UtilityMatrix + permute(repmat(beta.*Value*MarkovMatrix',[1,1,AGridN]),[1,3,2]), [], 1), AGridN, 9);
        else
            if mod(iter,10)==0
                [NewValue, Index] = max(UtilityMatrix + permute(repmat(beta.*Value*MarkovMatrix',[1,1,AGridN]),[1,3,2]), [], 1);
                NewValue = reshape(NewValue, AGridN, 9);
                Index = reshape(Index, AGridN, 9);
            else
                for i=1:AGridN
                    for j=1:ZGridN
                        NewValue(i,j) = UtilityMatrix(Index(i, j), i, j) + ...
                            beta.*Value(Index(i, j),:)*MarkovMatrix(j,:)';
                    end
                end 
            end
        end
        error = max(abs(Value - NewValue),[],'all',"linear");
        Value = NewValue;
        iter = iter+1;
    end
end