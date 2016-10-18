function getFRF(self)
%% getFRF via residues
% 
% 
% 
% author: john devitis
% create date: 18-Oct-2016 17:51:16
	
    if isempty(self.AA)
        self.getResidues();
    end
    
    % get FRF via residues
	self.HH = zeros(self.no,self.ni,self.ns);
	
    for ii = 1:self.ne              % loop modes
        for jj = 1:self.no          % loop outputs
            for kk = 1:self.ni      % loop inputs
                out = self.out(jj);
                in = self.in(kk);
                for ll = 1:self.ns  % loop spectral lines
                    % form [HH] - add mode ii contribution -> [no,ni,ns]
                    %           - add complex conjugate
                    tt = self.AA(out,in,ii)./(1j*self.w(ll)-self.root(ii)) ...
                        + conj(self.AA(out,in,ii))./(1j*self.w(ll)-conj(self.root(ii)));
                    % add each mode for total response
                    self.HH(jj,kk,ll) = self.HH(jj,kk,ll)+tt;
                end
            end
        end
    end
	
end
