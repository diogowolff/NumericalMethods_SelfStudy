function [V] = EigenInvariantDist(PolicyIndex, MarkovMatrix)
    TransitionMatrix = zeros(numel(PolicyIndex), numel(PolicyIndex));
    SizeAGrid = size(PolicyIndex, 1);

    for i=1:size(MarkovMatrix, 1)
        Aux = zeros(SizeAGrid, SizeAGrid);
        Aux(sub2ind([SizeAGrid, SizeAGrid], 1:SizeAGrid, PolicyIndex(:,1)')) = 1;
        TransitionMatrix(((i-1)*SizeAGrid+1):(i*SizeAGrid), ...
            :) = kron(MarkovMatrix(i,:), Aux);
    end
    
    [V] = eigs(TransitionMatrix, 1);
end

