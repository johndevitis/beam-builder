function [K, M] = beam_builder(E, I, L, ne, mbar, mchk)
%% function [K, M] = beambuilder(E, I, L, ne, mbar, mchk, condind)
%
% Assembles SS beam with consistent or lumped mass. Rotational DOF are
% removed (but actions kept) using static condensation (stiffness terms)
% and guyan reduction (mass terms). 
% Mass can be lumped or consistent (chk = 1 -> lumped). 
%
% condind = dof index for static condensation
% mchk = 1 -> lumped mass
%
% jdv 7/25/13; 08182015; 08192015; 08202015; 09142015

% setup
l = L/ne;  % element length
nn = ne+1; % number of nodes 

% get element matrices
[k,m] = getelement(E,I,l,mbar,mchk);

% assemble global matrices
[K,M] = assemble(k,m,nn);

%remove support dof
bc = [1 2 length(K)-1:length(K)]; 
K = removerows(K,'ind',bc);
K = removerows(K','ind',bc); % transpose
M = removerows(M,'ind',bc);
M = removerows(M','ind',bc); % transpose

end

% choose element k and m - translation and rotation
function [k,m] = getelement(E,I,l,mbar,mchk)
    % Element stiffness matrix [k] - bernoulli beam
    k = (E*I/l^3) * [  12,   6*l,  -12,   6*l ;...
                      6*l, 4*l^2, -6*l, 2*l^2 ;...
                      -12,  -6*l,   12,  -6*l ;...
                      6*l, 2*l^2, -6*l, 4*l^2 ];
    % Element Mass Matrix
    if mchk == 1 % check lumped for cont
        % element continuous mass matrix [m]
        m = (mbar*l/2) * [  1,0,0,0;...
                            0,0,0,0;...
                            0,0,1,0;...
                            0,0,0,0];         
    else
        % element continuous mass matrix [m]
        m = (mbar*l/420) * [  156,   22*l,    54,  -13*l ;...
                             22*l,  4*l^2,  13*l, -3*l^2 ;...
                               54,   13*l,   156,  -22*l ;...
                            -13*l, -3*l^2, -22*l,  4*l^2 ];        
    end
end


% Assemble K & M
function [K,M] = assemble(k,m,nn)
% k - element stiffness
% m - element mass
% nn - number of nodes
    % pre-allocate 
    K = zeros(2*nn,2*nn);
    M = zeros(size(K));
    % loop elements - (thanks, nate)
    for ii = 1:nn-1
        K((2*ii-1):(2*ii+2),(2*ii-1):(2*ii+2)) = K((2*ii-1):(2*ii+2),(2*ii-1):(2*ii+2)) + k;
        M((2*ii-1):(2*ii+2),(2*ii-1):(2*ii+2)) = M((2*ii-1):(2*ii+2),(2*ii-1):(2*ii+2)) + m;
    end
end






