function assemble(self)
%% assemble K & M
%
% k - element stiffness
% m - element mass
% nn - number of nodes
% 
% 
% author: john devitis
% create date: 18-Oct-2016 14:22:02
	

% pre-allocate 
K = zeros(2*self.nn,2*self.nn);
M = zeros(size(K));
% loop elements - (thanks, nate)
for ii = 1:self.nn-1
    K((2*ii-1):(2*ii+2),(2*ii-1):(2*ii+2)) = K((2*ii-1):(2*ii+2),...
                                                (2*ii-1):(2*ii+2)) + ...
                                                self.k;
    M((2*ii-1):(2*ii+2),(2*ii-1):(2*ii+2)) = M((2*ii-1):(2*ii+2),...
                                                (2*ii-1):(2*ii+2)) + ...
                                                self.m;
end

% assign to self
self.K = K;
self.M = M;

end
