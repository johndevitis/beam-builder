function [Hs,hh] = vibsFRF(AA,root,in,out,w)
%% Vibs Example - Construct FRF - Mode Contributions
%
% logmag plot of Hpq with dashed mode contributions
% complex conjugate solution
%
% jdv 06232015; 07232015; 08162015

%% Setup

% spatial sampling
no = length(out);
ni = length(in);
ne = length(root);

%% Get FRF

% get individual contributions
hh = zeros(length(w),ne);
for ii = 1:ne
    for jj = 1:length(w)
        % mode ii contribution to FRF
        %   note: add complex conjugate
        hh(jj,ii) = AA(out,in,ii)./(1j*w(jj)-root(ii)) + ...
                    conj(AA(out,in,ii))./(1j*w(jj)-conj(root(ii))); 
    end
end

% sum for full response
Hs = sum(hh,2);

