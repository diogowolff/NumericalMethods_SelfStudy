function [mu] = EigenInvariantDist(PolicyIndex, MarkovMatrix)
    TransitionMatrix = zeros(numel(PolicyIndex), numel(PolicyIndex));
    SizeAGrid = size(PolicyIndex, 1);

    for i=1:size(MarkovMatrix, 1)
        Aux = zeros(SizeAGrid, SizeAGrid);
        Aux(sub2ind([SizeAGrid, SizeAGrid], 1:SizeAGrid, PolicyIndex(:,1)')) = 1;
        TransitionMatrix(((i-1)*SizeAGrid+1):(i*SizeAGrid), ...
            :) = kron(MarkovMatrix(i,:), Aux);
    end
    
    mu = (eye(size(TransitionMatrix, 1)) - TransitionMatrix + ones(size(TransitionMatrix)))\ones(size(TransitionMatrix, 1),1);
end

