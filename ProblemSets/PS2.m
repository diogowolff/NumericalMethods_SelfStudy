%% Question 3

clear
addpath('../AuxiliaryMethods');

beta = .987;
zeta = 1/beta - 1;
mu = 2;
alpha = 1/3;
delta = .012;
rho = .95;
sigma = .007;

% Generate discretized uncertainty

[ShockGrid, MarkovMatrix] = TauchenDiscretizer(7, 3, 0, rho, sigma);

SteadyState = ((zeta-delta)/alpha)^(1/(alpha-1));
StockGrid = linspace(.75*SteadyState, 1.25*SteadyState, 500);
[KGrid, KNewGrid, TFPGrid] = meshgrid(StockGrid, StockGrid, exp(ShockGrid));

UtilityMatrix = ((TFPGrid.*KGrid.^alpha + (1-delta).*KGrid - KNewGrid).^(1-mu)-1)./(1-mu);

Value = zeros(500, 7);
NewValue = Value;
error = 100;
tol = 10e-6;
iter = 0;

tic
while error > tol
    for i=1:500
        for j=1:7
            PossibleNewValues = UtilityMatrix(i, i:500, j) + ...
                beta.*Value(i:500, :)*MarkovMatrix(j,:)';
            NewValue(i,j) = max(PossibleNewValues,[],'all');
        end
    end
    error = max(abs(Value - NewValue),[],'all');
    Value = NewValue;
    iter = iter+1;
end
toc

plot(Value(:,1))