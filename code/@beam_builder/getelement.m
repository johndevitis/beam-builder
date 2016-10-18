function getelement(self,bm)
%% getelement
% 
% 
% 
% author: john devitis
% create date: 18-Oct-2016 14:19:50
	
% choose element k and m - translation and rotation
    l = self.l; 
    E = bm.E;
    I = bm.I;
    
    % Element stiffness matrix [k] - bernoulli beam
    self.k = (E*I/l^3) * [  12,   6*l,  -12,   6*l ;...
                           6*l, 4*l^2, -6*l, 2*l^2 ;...
                           -12,  -6*l,   12,  -6*l ;...
                           6*l, 2*l^2, -6*l, 4*l^2 ];
    % Element Mass Matrix
    if self.mchk == 1 % check lumped for cont
        % element continuous mass matrix [m]
        self.m = (bm.mbar*l/2) * [  1,0,0,0;...
                                      0,0,0,0;...
                                      0,0,1,0;...
                                      0,0,0,0];         
    else
        % element continuous mass matrix [m]
        self.m = (bm.mbar*l/420) * [  156,   22*l,    54,  -13*l ;...
                                       22*l,  4*l^2,  13*l, -3*l^2 ;...
                                         54,   13*l,   156,  -22*l ;...
                                      -13*l, -3*l^2, -22*l,  4*l^2 ];        
    end
end
	

