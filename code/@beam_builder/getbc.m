function getbc(self,macro)
%% getbc
% 
% get boundary condition macro
% 
% author: john devitis
% create date: 18-Oct-2016 15:49:15
	
    % error screen null entry
	if nargin < 2; macro = 'simple'; end;
    
	if ischar(macro)
        switch macro
            case 'simple' 
                self.bc = [1 2 length(self.K)-1:length(self.K)];
            case 'ss'
                self.bc = [1 2 length(self.K)-1:length(self.K)];
        end
    end
	
	
end
