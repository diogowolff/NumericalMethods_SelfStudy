function [ShockGrid, MarkovMatrix] = TauchenDiscretizer(n,scale,mu,rho,sigma)
    UpperBound = scale*sigma/sqrt(1-rho^2);
    LowerBound = -UpperBound;

    ShockGrid = linspace(LowerBound, UpperBound, n);
    dTheta = ShockGrid(2)-ShockGrid(1);


    MarkovMatrix = zeros(n, n);

    MarkovMatrix(:,1) = normcdf(...
        (LowerBound + dTheta/2 - (1-rho)*mu - rho.*ShockGrid)./sigma);

    MarkovMatrix(:,n) = 1 - normcdf(...
        (UpperBound - dTheta/2 - (1-rho)*mu - rho.*ShockGrid)./sigma);

    for j = 2:(n-1)
        MarkovMatrix(:,j) = normcdf(...
            (ShockGrid(j) + dTheta/2 - (1-rho)*mu - rho.*ShockGrid)./sigma) - ...
        normcdf(...
            (ShockGrid(j) - dTheta/2 - (1-rho)*mu - rho.*ShockGrid)./sigma);
    end
end

