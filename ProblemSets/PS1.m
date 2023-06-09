clear

addpath('../AuxiliaryMethods');

%% Question 1

rho = .95;
sigma = .007;
mu = 0;
n = 9;
scale = 3; % scale as in slides

[ShockGrid1, MarkovMatrix1] = TauchenDiscretizer(n, scale, mu, rho, sigma);

%% Question 2

[ShockGrid2, MarkovMatrix2] = RouwenhorstDiscretizer(n, rho, sigma);

%% Question 3

errorDraws = normrnd(0, .007, 10000, 1);
simulatedAR = zeros(10000,1);

simulatedAR(1)=errorDraws(1);

for i=2:10000
    simulatedAR(i) = rho*simulatedAR(i-1) + errorDraws(i);
end


tauchenStateSimulation = simulate(dtmc(MarkovMatrix1), 10000);
tauchenSimulation = ShockGrid1(tauchenStateSimulation);

rouwenhorstStateSimulation = simulate(dtmc(MarkovMatrix2), 10000);
rouwenhorstSimulation = ShockGrid2(rouwenhorstStateSimulation);

tiledlayout(2,1);

nexttile
plot(simulatedAR);
hold on
plot(tauchenSimulation);
title('Tauchen discretization of AR process');

nexttile
plot(simulatedAR);
hold on
plot(rouwenhorstSimulation);
title('Rouwenhorst discretization of AR process');

%The discrete models seem to work fine.

%% Question 4

tauchenEstimate = ar(tauchenSimulation, 1);
rouwenhorstEstimate = ar(rouwenhorstSimulation, 1);

%They're both pretty good.

%% Question 5

rho = .99;

[ShockGrid3, MarkovMatrix3] = TauchenDiscretizer(n, scale, mu, rho, sigma);
[ShockGrid4, MarkovMatrix4] = RouwenhorstDiscretizer(n, rho, sigma);

simulatedAR2(1)=errorDraws(1);

for i=2:10000
    simulatedAR2(i) = rho*simulatedAR2(i-1) + errorDraws(i);
end

tauchenStateSimulation2 = simulate(dtmc(MarkovMatrix3), 10000);
tauchenSimulation2 = ShockGrid1(tauchenStateSimulation2);

rouwenhorstStateSimulation2 = simulate(dtmc(MarkovMatrix4), 10000);
rouwenhorstSimulation2 = ShockGrid2(rouwenhorstStateSimulation2);

tiledlayout(2,1);

nexttile
plot(simulatedAR2);
hold on
plot(tauchenSimulation2);

nexttile
plot(simulatedAR2);
hold on
plot(rouwenhorstSimulation2);

%Huh, as expected, Tauchen's method seems to break down with a rho close to
%1. Clearly, it fails to reach the peaks and valleys that the true process
%attains. Rouwenhorst's seems to still work.

tauchenEstimate2 = ar(tauchenSimulation2, 1);
rouwenhorstEstimate2 = ar(rouwenhorstSimulation2, 1);

%Tauchen overshoots a little but works well enough; Rouwenhorst is more
%accurate. Predicting the rho seems to be fine in both cases, just the
%actual simulations that are starkly different.