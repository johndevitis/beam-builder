function [hs,h] = vibsIRF(AA,root,in,out,fs,l)
%% Vibs Example - Impulse Resopnse Function
%
%
% jdv 06232015; 07232015; 08162015

%% setup

% time sampling parameters
dt = 1/fs;
t  = 0:dt:l-dt;

%% get impulse response function

% get individual contributions
ne = length(root);
h = zeros(length(t),ne);
for ii = 1:ne
    for jj = 1:length(t)
        % mode ii contribution to IRF
        h(jj,ii) = AA(out,in,ii)*exp(root(ii)*t(jj)) + ...
                    conj(AA(out,in,ii))*exp(conj(root(ii))*t(jj)); 
    end
end

% sum for full response
hs = sum(h,2);