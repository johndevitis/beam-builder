%% Publish Build Details and Analysis
%
% jdv 11132015

%% HTML Options

% root write path
pname = 'C:\Users\John\Dropbox\Thesis Docs - Dropbox\Numerical Models\5DOF Beam\Publishing';
mkdir(pname);
h = htm(pname);  % create instance
h.start;         % open file

% text
h.h1('This is a test Header1')
h.h2('This is a test Header2')
h.comment('this is a comment')
h.dev('this is a dev')






