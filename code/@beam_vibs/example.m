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
    bv.K = bb.K;
    bv.M = bb.M;
    bv.eig();

    bv.dampr = ones(size(bv.W))*.05; % 5% critical damping
    
    bv.ns = 2^9; 
    bv.freqbnds = [-150 150];
    
    bv.in = 1:5;
    bv.out = 1:5;
    
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
