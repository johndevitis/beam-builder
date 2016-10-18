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
	
	
end
