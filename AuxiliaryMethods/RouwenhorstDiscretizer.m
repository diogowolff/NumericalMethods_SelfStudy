function [ShockGrid,MarkovMatrix] = RouwenhorstDiscretizer(n,rho,sigma)
    UpperBound = sqrt(sigma^2/(1-rho^2))*sqrt(n-1);
    LowerBound = -UpperBound;

    ShockGrid = linspace(LowerBound, UpperBound, n);
    
    MarkovMatrix = zeros(n, n);
    
    p = (1+rho)/2;
    
    MarkovMatrix(1:2, 1:2) = [p, 1-p; 1-p, p];
    
    for i=3:n
       zero = zeros(i-1, 1);
       MarkovMatrix(1:i, 1:i) = p*[MarkovMatrix(1:(i-1),1:(i-1)), zero; ...
           zero', 0] + (1-p)*[zero, MarkovMatrix(1:(i-1),1:(i-1)); ...
           0, zero'] + (1-p)*[zero', 0; MarkovMatrix(1:(i-1),1:(i-1)), zero] + ...
           p*[0, zero'; zero, MarkovMatrix(1:(i-1),1:(i-1))];
    end
    
    MarkovMatrix = normr(MarkovMatrix);
end

