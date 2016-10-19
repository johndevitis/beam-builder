function example()
%% example
% 
% 
% 
% author: john devitis
% create date: 18-Oct-2016 16:49:31
	
    % use defaults from 5dof example in bb
	bb = beam_builder.example5();
    
	% create instance of beam_vibs class and assign parameters
    bv = beam_vibs();
    
    % assign mass and stiffness
    bv.K = bb.K;
    bv.M = bb.M;
    
    % eigen-solution and pole sorting
    bv.eig();

    % assign percent critical damping
    bv.dampr = ones(size(bv.W))*.05; % 5% critical damping
    
    % assign spectral info
    bv.ns = 2^9;                % number of freq samples
    bv.freqbnds = [-150 150];   % -/+ freq bounds
    
    % frf spatial sampling1
    bv.in = 1:5;
    bv.out = 1:5;
    
    % get frf
    bv.getFRF();
    
    
    %% H11 - FRF
    in = 1;
    out = 1;
    [Hs,hh] = bv.vibsFRF(bv.AA,bv.root,in,out,bv.w);
    bv.vibsFRFplot(Hs,hh,in,out,bv.w);
    
    %% H11 - Phase
    bv.vibsPhaseplot(Hs,hh,in,out,bv.w);


    %% H11 - IRF
    fs = 200; % sampling freq
    l = .5;    % length [sec] 
    [hs,h] = bv.vibsIRF(bv.AA,bv.root,in,out,fs,l);
    bv.vibsIRFplot(hs,h,in,out,fs,l);
    
    
    %% plot mode shapes
    fh = bb.plot();
    bv.plotShape(bb)
    
end
