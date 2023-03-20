clear

addpath('../AuxiliaryMethods');

%% item a)

rho = .9;
sigma = .01;

[ShockGrid, MarkovMatrix] = TauchenDiscretizer(9, 3, 0, rho, sigma);

%% item b)

[Value, UtilityMatrix, MarkovMatrix, AssetGrid] = ExperimentalHuggett([.96, 1.0001, .9, .01], 999*(1/beta - 1)/1000, 1000, 9);

%% item c)

[test, test2] = PolicyHuggett(.96, Value, AssetGrid, 9, UtilityMatrix, MarkovMatrix);
test4 = InvariantDistribution(test2, MarkovMatrix);
test5 = EigenInvariantDist(test2, MarkovMatrix);

space = linspace(0, 1/beta-1, 100);
[test7(1), dist] = AssetDemandFunction(0, Value);
for i=2:100
    [test7(i), dist] = AssetDemandFunction(space(i), Value, dist);
    disp(i)
end

test7 = HuggettEquilibrium([.96, 1.0001, .9, .01], 1000, 9);