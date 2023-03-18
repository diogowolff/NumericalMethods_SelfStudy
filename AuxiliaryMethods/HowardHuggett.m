function [Value, UtilityMatrix, MarkovMatrix, AssetGrid] = HowardHuggett(param, r, AGridN, ZGridN, StartingGuess)
    
    beta = param(1);
    gamma = param(2);
    rho = param(3);
    sigma = param(4);
    
    

    [ShockGrid, MarkovMatrix] = TauchenDiscretizer(ZGridN, 3, 0, rho, sigma);
    
    phi = (exp(ShockGrid(1))/(1-beta)); 

    AssetGrid = linspace(-phi, phi, AGridN);
    [AGrid, ANewGrid, TFPGrid] = meshgrid(AssetGrid, AssetGrid, exp(ShockGrid));

    UtilityMatrix = ((exp(TFPGrid) + (1+r).*AGrid - ANewGrid).^(1-gamma)-1)./(1-gamma);
    UtilityMatrix(exp(TFPGrid) + (1+r).*AGrid - ANewGrid <=0 ) = -inf;

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
            for i=1:AGridN
                for j=1:ZGridN
                    PossibleNewValues = UtilityMatrix(:, i, j) + ...
                        beta.*Value*MarkovMatrix(j,:)';
                    NewValue(i,j) = max(PossibleNewValues,[],'all',"linear");
                end
            end
        else
            if mod(iter,10)==0
                for i=1:AGridN
                    for j=1:ZGridN
                        PossibleNewValues = UtilityMatrix(:, i, j) + ...
                            beta.*Value*MarkovMatrix(j,:)';
                        [NewValue(i,j), Index(i, j)]= max(PossibleNewValues,[],'all',"linear");
                    end
                end     
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

