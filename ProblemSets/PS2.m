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

SteadyState = ((zeta+delta)/alpha)^(1/(alpha-1));
StockGrid = linspace(.75*SteadyState, 1.25*SteadyState, 500);
[KGrid, KNewGrid, TFPGrid] = meshgrid(StockGrid, StockGrid, exp(ShockGrid));

UtilityMatrix = ((TFPGrid.*KGrid.^alpha + (1-delta).*KGrid - KNewGrid).^(1-mu)-1)./(1-mu);
UtilityMatrix(TFPGrid.*KGrid.^alpha + (1-delta).*KGrid - KNewGrid <=0 ) = -inf;

%% Brute force

%Starting with u(SS)/(1-beta) as it is a good guess
Value = zeros(500, 7) + ((SteadyState^alpha - delta*SteadyState)^(1-mu)-1)/((1-mu)*(1-beta));
NewValue = Value;
error = 100;
tol = 10e-6;
iter = 0;

tic
while error > tol
    for i=1:500
        for j=1:7
            PossibleNewValues = UtilityMatrix(:, i, j) + ...
                beta.*Value*MarkovMatrix(j,:)';
            NewValue(i,j) = max(PossibleNewValues,[],'all');
        end
    end
    error = max(abs(Value - NewValue),[],'all');
    Value = NewValue;
    iter = iter+1;
end
toc

%7.15 seconds

%% Using monotonicity

%Starting with u(SS)/(1-beta) as it is a good guess
Value = zeros(500, 7) + ((SteadyState^alpha - delta*SteadyState)^(1-mu)-1)/((1-mu)*(1-beta));
NewValue = Value;
error = 100;
tol = 10e-6;
iter = 0;
ChoiceGivenState = zeros(7,1);


tic
while error > tol
    for i=1:500
        if i==1
            for j=1:7
            PossibleNewValues = UtilityMatrix(:, i, j) + ...
                beta.*Value*MarkovMatrix(j,:)';
            [NewValue(i,j), ChoiceGivenState(j)] = max(PossibleNewValues,[],'all');
            end
        else   
            for j=1:7
            PossibleNewValues = UtilityMatrix(ChoiceGivenState(j):500, i, j) + ...
                beta.*Value(ChoiceGivenState(j):500, :)*MarkovMatrix(j,:)';
            [NewValue(i,j), ChoiceGivenState(j)] = max(PossibleNewValues,[],'all');
            end
        end
    end
    error = max(abs(Value - NewValue),[],'all');
    Value = NewValue;
    iter = iter+1;
end
toc

% 12.6 seconds; possibly would be faster if I were to parallelize the code.

%Some thoughts: the evaluation over each point of the shock grid is
%essentially independent, being connected only through the markov matrix.
%If I reverse the order of iteration to shock(stock) I can just put a
%parfor on the outside loop and it'll all be parallel.

%% Policy function

Index = zeros(500, 7);
Choice = zeros(500, 7);

for i=1:500
    for j=1:7
        PossibleValues = UtilityMatrix(:, i, j) + ...
                beta.*Value*MarkovMatrix(j,:)';
        [~, Index(i, j)] = max(PossibleValues,[],'all');
        Choice(i,j) = StockGrid(Index(i, j));
    end
end

plot(Choice(220:280,:));

% Doing the policy like this was a really bad idea. Creating an actual
% function would be much better, might do it later.


%% Euler errors

ConsumptionToday = (reshape(TFPGrid(1,:,:).*KGrid(1,:,:).^alpha + (1-delta).*KGrid(1,:,:),[500,7]) - Choice);

Index2 = zeros(500, 7);
Choice2 = zeros(500, 7);
for i=1:500
    for j=1:7
        Index2(i, j) = Index(Index(i,j), j);
        Choice2(i,j) = StockGrid(Index2(i, j));
    end
end

MarginalUtilityTomorrow = (reshape(TFPGrid(1,:,:),[500,7]).*Choice.^alpha + (1-delta).*Choice - Choice2).^(-mu).*...
    (alpha.*(reshape(TFPGrid(1,:,:),[500,7]).*Choice.^(alpha-1)) + (1-delta));
ExpecMgUtTomorrow = MarginalUtilityTomorrow*MarkovMatrix';
EEE = log10(abs(1 - (beta.*ExpecMgUtTomorrow).^(-1/mu)./ConsumptionToday));

plot(EEE);

%It works!



%% Question 4

%Starting with u(SS)/(1-beta) as it is a good guess
Value = zeros(500, 7) + ((SteadyState^alpha - delta*SteadyState)^(1-mu)-1)/((1-mu)*(1-beta));
NewValue = Value;
error = 100;
tol = 10e-6;
iter = 0;



tic
while error > tol
    if iter<100
        for i=1:500
            for j=1:7
                PossibleNewValues = UtilityMatrix(:, i, j) + ...
                    beta.*Value*MarkovMatrix(j,:)';
                NewValue(i,j) = max(PossibleNewValues,[],'all');
            end
        end
    else
        if mod(iter,10)==0
            for i=1:500
                for j=1:7
                    PossibleNewValues = UtilityMatrix(:, i, j) + ...
                        beta.*Value*MarkovMatrix(j,:)';
                    [NewValue(i,j), Index(i, j)]= max(PossibleNewValues,[],'all');
                end
            end     
        else
            for i=1:500
                for j=1:7
                    NewValue(i,j) = UtilityMatrix(Index(i, j), i, j) + beta.*Value(Index(i, j),:)*MarkovMatrix(j,:)';
                end
            end 
        end
    end
    display(error)
    error = max(abs(Value - NewValue),[],'all');
    Value = NewValue;
    iter = iter+1;
end
toc