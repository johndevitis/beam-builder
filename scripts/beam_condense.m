 %% static condensation and guyan reduction to remove dof
function [Kr,Mr] = beam_condense(K,M,rm,mchk)
% function Kr = condense(K,M,rm,mchk)
% K  - stiffness matrix
% M  - mass matrix
% rm - dof index to remove
% mchk - 1 for lumped mass - preserves diagonal quantities
% mchk - 0 for lumped or consistent 
%
% jdv 09142015; 09152015
    
    % error screen dof index to remove
    if size(rm,1) > size(rm,2)
        rm = rm';             
    end
    
    % total dof
    nn = length(K);           
    
    % find index to keep
    aa = setdiff(1:nn,rm);
    
    % re-order global matrices
    KK = K([aa rm],[aa rm]);
    
    % get index of portion to be removed
    ll = 1:length(KK)-length(rm);   % left to keep
    rr = (length(aa)+1):length(KK); % right to remove
    
    % get partitions of K matrix
    krr = KK(ll,ll); % retained
    kcr = KK(rr,ll); % coupled
    kcc = KK(rr,rr); % condensed
    
    % get dof transform matrix
    dd = -(kcc^-1)*kcr;   % denom of transform
    ii = eye(size(krr));  % numer of transform 
    T = [ii;dd];          % transform matrix
    
    % static condensation to remove dof
    Kr = T'*KK*T;  
    
    % check for mass option
    if ~isempty(M) || ~isempty(mchk)
        % check mass type
        if mchk == 1 
            % lumped mass - just remove zeros
            M = removerows(M,'ind',rm);  
            Mr = removerows(M','ind',rm); % transpose
        else 
            % guyan reduction - continuous or lumped (smears mass)
            MM = M([aa rm],[aa rm]);    % re-order global matrix
            Mr = T'*MM*T;               % remove dof   
        end
    else
        Mr = [];
    end
    
end