function [Ur,Sr,Vr,rank] = truncSVD(X,rank,tSize)
    [nObs, nStates] = size(X);
    
    [U,S,V] = svd(X,'econ');
    
    if isempty(rank) || rank <= 0
        beta = tSize / nObs;
        sigma = diag(S);
        tau = optimal_SVHT_coef(beta,0) * median(sigma);

        rank = find(sigma > tau,1,'last');
    elseif rank >= nStates
       rank = nStates;
    end
    
    Ur = U(:,1:rank);
    Sr = S(1:rank,1:rank);
    Vr = V(:,1:rank);
end

