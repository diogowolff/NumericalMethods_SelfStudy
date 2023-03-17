function [Value] = HowardVFI(param,KGridN,ZGridN, StartingGuess)
    beta = param(1);
    zeta = param(2);
    mu = param(3);
    alpha = param(4);
    delta = param(5);
    rho = param(6);
    sigma = param(7);
    
    [ShockGrid, MarkovMatrix] = TauchenDiscretizer(ZGridN, 3, 0, rho, sigma);
    
    SteadyState = ((zeta+delta)/alpha)^(1/(alpha-1));
    StockGrid = linspace(.75*SteadyState, 1.25*SteadyState, KGridN);
    [KGrid, KNewGrid, TFPGrid] = meshgrid(StockGrid, StockGrid, exp(ShockGrid));

    UtilityMatrix = ((TFPGrid.*KGrid.^alpha + (1-delta).*KGrid - KNewGrid).^(1-mu)-1)./(1-mu);
    UtilityMatrix(TFPGrid.*KGrid.^alpha + (1-delta).*KGrid - KNewGrid <=0 ) = -inf;
    
    
    % Starting guess; if not given then start with K_ss
    
    if ~exist('StartingGuess','var')
         %Starting with u(SS)/(1-beta) as it is a good guess
        Value = zeros(KGridN, ZGridN) + ...
            ((SteadyState^alpha - delta*SteadyState)^(1-mu)-1)/((1-mu)*(1-beta));
    else
        Value = StartingGuess;
    end
    
    
    
    NewValue = Value;
    error = 100;
    tol = 10e-6;
    iter = 0;



    
    while error > tol
        if iter<50
            for i=1:KGridN
                for j=1:ZGridN
                    PossibleNewValues = UtilityMatrix(:, i, j) + ...
                        beta.*Value*MarkovMatrix(j,:)';
                    NewValue(i,j) = max(PossibleNewValues,[],'all',"linear");
                end
            end
        else
            if mod(iter,10)==0
                for i=1:KGridN
                    for j=1:ZGridN
                        PossibleNewValues = UtilityMatrix(:, i, j) + ...
                            beta.*Value*MarkovMatrix(j,:)';
                        [NewValue(i,j), Index(i, j)]= max(PossibleNewValues,[],'all',"linear");
                    end
                end     
            else
                for i=1:KGridN
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

