function getResidues(self)
%% getResidues
% 
% 
% 
% author: john devitis
% create date: 18-Oct-2016 17:38:26
	AA = zeros(self.no,self.ni,self.ne);
	for ii = 1:self.ne % loop modes
        self.AA(:,:,ii) = self.Qr(ii) * self.V(:,ii) * self.V(:,ii)';
    end	
end
