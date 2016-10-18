function bb = example5()
%% example
% 
% bh - beam_builder object handle
% 
% author: john devitis
% create date: 18-Oct-2016 14:30:28
%   * original edits: jdv 06232015; 07232015; 08162015; 11122015
	
	fprintf('Vibs Example - 5DOF Lumped Mass Shear Beam\n');
    
    %% Create instance of beam section object
    bm = beam_rect() % use defaults
    
    
    %% Create instance of beam builder object and input build parameters
    bb = beam_builder();
    
    %%
    % build parameters:
    %   nodes: 1-2-3-4-5-6-7 
    %   dof:  bc-1-2-3-4-5-bc
    bb.L = 100*12; % total length 100ft -> in
    bb.nn = 7; % total nodes
    
    % get element k, m matrices for beam section, bm
    bb.getelement(bm); 
    
    % assemble global K, M matrices
    bb.assemble();
    
    % use getbc macro method for getting bc indices
    bb.getbc('simple'); % this populates the bc property of beam_builder

    % apply boundary conditions
    bb.applybc(); 
    
    
    %% Make shear beam - remove rotation dof
	ind = 2:2:length(bb.K);
    bb.applybc(ind);
    
    
    bb % display contents  
    
	% plot beam
    bb.plot();    
    
end













