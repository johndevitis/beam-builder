function eig(self)
%% eig
% 
% 
% 
% author: john devitis
% create date: 18-Oct-2016 16:49:18
	
    [V,D] = eig(self.K,self.M); % solve
    [val,ind] = sort(diag(D));  % sort eigenvalues
    self.V = V(:,ind);          % apply sort to eigenvectors
    self.W = sqrt(val);         % [rad/sec]
    self.F = self.W/2/pi;       % [hz]

    % decouple system matrices
    self.Mr = self.V'*self.M*self.V; % modal mass, Mr
    self.Kr = self.V'*self.K*self.V; % modal stiffness, Kr

    % form mass normalized modeshapes
    self.Vn = zeros(size(self.V));
%     self.ne = size(self.V,2);     % number of effective modes
    for ii = 1:self.ne
        self.Vn(:,ii) = self.V(:,ii)/sqrt(self.Mr(ii,ii));
    end
	
end
