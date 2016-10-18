function example()
%% example
% 
% 
% 
% author: john devitis
% create date: 18-Oct-2016 16:49:31
	
    % use defaults from 5dof example in bb
	bb = beam_builder.example5();
    
	% create instance of beam_vibs class
    bv = beam_vibs();
    bv.K = bb.K;
    bv.M = bb.M;
    bv.eig();
	
    bv.dampr = ones(size(bv.W))*.05; % 5% critical damping
    
    bv.ns = 2^9; 
    bv.freqbnds = [-150 150];
    
    bv.in = 1:5;
    bv.out = 1:5;
    
    bv.getFRF();
    
	
end
