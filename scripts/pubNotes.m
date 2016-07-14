%% Publishing 5DOF Beam
%
% jdv 11152015



%% Publish Build Details and Analysis

% path to publish to
pname = 'C:\Users\John\Dropbox\Thesis Docs - Dropbox\Numerical Models\5DOF Beam\Publishing\Design';
mkdir(pname);
% publishing options 
opts = struct('format','html',...
              'outputDir',pname);          
% publish code
publish('beam_5DOF_design.m',opts);
close all


pname = 'C:\Users\John\Dropbox\Thesis Docs - Dropbox\Numerical Models\5DOF Beam\Publishing\Design';
% publishing options 
opts = struct('format','pdf',...
              'outputDir',pname);          
% publish code
publish('beam_5DOF_design.m',opts);
close all

pname = 'C:\Users\John\Dropbox\Thesis Docs - Dropbox\Numerical Models\5DOF Beam\Publishing\Design';
% publishing options 
opts = struct('format','doc',...
              'outputDir',pname);          
% publish code
publish('beam_5DOF_design.m',opts);
close all


%% Publish THMPR - Scenario 1


% path to publish to
pname = 'C:\Users\John\Dropbox\Thesis Docs - Dropbox\Numerical Models\5DOF Beam\Publishing\Scenario 1';
mkdir(pname);
% publishing options - HTML
opts = struct('format','html',...
              'outputDir',pname);          
% publish code
publish('beam_5DOF_thmpr_1.m',opts);
close all


pname = 'C:\Users\John\Dropbox\Thesis Docs - Dropbox\Numerical Models\5DOF Beam\Publishing\Scenario 1';
% publishing options
opts = struct('format','doc',...
              'outputDir',pname);          
% publish code
publish('beam_5DOF_thmpr_1.m',opts);
close all




%% Publish THMPR - Scenario 2


% path to publish to
pname = 'C:\Users\John\Dropbox\Thesis Docs - Dropbox\Numerical Models\5DOF Beam\Publishing\Scenario 2';
mkdir(pname);
% publishing options - HTML
opts = struct('format','html',...
              'outputDir',pname);          
% publish code
publish('beam_5DOF_thmpr_2.m',opts);
close all


pname = 'C:\Users\John\Dropbox\Thesis Docs - Dropbox\Numerical Models\5DOF Beam\Publishing\Scenario 2';
% publishing options - PDF 
opts = struct('format','pdf',...
              'outputDir',pname);          
% publish code
publish('beam_5DOF_thmpr_2.m',opts);
close all

pname = 'C:\Users\John\Dropbox\Thesis Docs - Dropbox\Numerical Models\5DOF Beam\Publishing\Scenario 2\Doc';
% publishing options - PDF 
opts = struct('format','doc',...
              'outputDir',pname);          
% publish code
publish('beam_5DOF_thmpr_2.m',opts);
close all



%% Publish THMPR - Scenario 3

% path to publish to
pname = 'C:\Users\John\Dropbox\Thesis Docs - Dropbox\Numerical Models\5DOF Beam\Publishing\Scenario 3';
mkdir(pname);
% publishing options 
opts = struct('format','html',...
              'outputDir',pname);          
% publish code
publish('beam_5DOF_thmpr_3.m',opts);
close all


pname = 'C:\Users\John\Dropbox\Thesis Docs - Dropbox\Numerical Models\5DOF Beam\Publishing\Scenario 3';
% publishing options 
opts = struct('format','pdf',...
              'outputDir',pname);          
% publish code
publish('beam_5DOF_thmpr_3.m',opts);
close all


pname = 'C:\Users\John\Dropbox\Thesis Docs - Dropbox\Numerical Models\5DOF Beam\Publishing\Scenario 3';
% publishing options 
opts = struct('format','doc',...
              'outputDir',pname);          
% publish code
publish('beam_5DOF_thmpr_3.m',opts);
close all



%% Publish THMPR - Scenario full

% path to publish to
pname = 'C:\Users\John\Dropbox\Thesis Docs - Dropbox\Numerical Models\5DOF Beam\Publishing\Scenario 4';
mkdir(pname);
% publishing options 
opts = struct('format','html',...
              'outputDir',pname);          
% publish code
publish('beam_5DOF_thmpr_3.m',opts);
close all


pname = 'C:\Users\John\Dropbox\Thesis Docs - Dropbox\Numerical Models\5DOF Beam\Publishing\Scenario 4';
% publishing options 
opts = struct('format','pdf',...
              'outputDir',pname);          
% publish code
publish('beam_5DOF_thmpr_3.m',opts);
close all


pname = 'C:\Users\John\Dropbox\Thesis Docs - Dropbox\Numerical Models\5DOF Beam\Publishing\Scenario 4';
% publishing options 
opts = struct('format','doc',...
              'outputDir',pname);          
% publish code
publish('beam_5DOF_thmpr_full.m',opts);
close all





