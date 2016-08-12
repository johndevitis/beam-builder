% remove rotational dof
function [K,M] = beam_removerotation(K,M,mchk)
% jdv 09152015

% check for mass type
if mchk == 1 % lumped mass
    % find zeros on diagonal
    rm = find(diag(M) == 0); 
    % if zeros - remove dof
    if ~isempty(rm)
        % condense any rotation dof
        [K,M] = beam_condense(K,M,rm,mchk);
    end
else % continuous
    % remove even index
    rm = 2:2:length(K);
    % condense rotations in stiffness matrix
    [K,M] = beam_condense(K,M,rm,mchk);
end

