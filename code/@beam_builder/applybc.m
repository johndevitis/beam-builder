function applybc(self,bc)
%% applybc
% 
% applies boundary condition(s) to beam
% 
% author: john devitis
% create date: 18-Oct-2016 15:19:18
    
    if nargin < 2
        bc = self.bc;
    end
    
	% remove support dof
    self.K = removerows(self.K,'ind',bc);
    self.K = removerows(self.K','ind',bc); % transpose
    self.M = removerows(self.M,'ind',bc); 
    self.M = removerows(self.M','ind',bc); % transpose
	
end
