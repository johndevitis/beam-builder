%% Solution of a proportionally damped (i.e. real normal modes) 2DOF beam 
%
% jdv 3/11/15; 3/12/15; 7/14/2015; 07/23/2015

%% Setup Structure 
% 2DOF bernoulli beam with discrete translational springs

% Design beam
b = 1;              % in
h = 12;             % in
I = b*h^3/12;       % in^4
E = 29e6;           % psi 
l = 100*12;         % ft -> in
grav = 386.4;       % in/sec^2;
mbar = b*h*.29;     % vol * density -> lbf
mbar = mbar/grav;   % lbm

% define discrete springs
k1 = E*I/l^3 * 300; % prev 300
k2 = E*I/l^3 * 300;

% form [M] [K]
M = mbar*l/2 * [ 1  0 ;  0 1];
k = E*I/l^3  * [ 12 -12; -12 12];

% Add springs to [K]
K = [k(1,1) + k1, k(1,2)    ;...
         k(2,1) , k(2,2)+k2];
     

%% Eigen analysis

% find weighted orthogonal vectors (orthonormal eigenvectors)
% and eigenvalues
[V,D] = eig(K,M); 

% get vector of un-dampted natural frequency
W = sqrt(diag(D)); % [rad/sec]

% decouple system 
Mr = V'*M*V; % modal mass, Mr - note: unitary mass
Kr = V'*K*V; % modal stiffness, Kr

% add proportional damping
dampr = .05;                 % damping ratio [5% critical damping]
dampf = -dampr*W;            % damping factor [rad/sec]
Wn = sqrt(W.^2 + dampf.^2);  % damped natural frequency [rad/sec]
root = dampf + 1j*Wn;        % form positive poles 


%% form transfer function from partial fraction expansion form

% define measurement temporal and spatial information
ns = 128;                   % number of spectral lines
no = 2;                     % number of outputs
ni = 2;                     % number of inputs
ne = 2;                     % number of effective modes
w = linspace(0,25,ns);      % [rad/sec]

% form residue matrices
Q1 = 1/(2j*root(1));    % mode 1 scaling for unit modal mass
Q2 = 1/(2j*root(2));    % mode 2 scaling for unit modal mass
A1 = Q1*V(:,1)*V(:,1)'; % [A] mode 1 residue matrix [no x ni] 
A2 = Q2*V(:,2)*V(:,2)'; % [A] mode 2 residue matrix [no x ni]

% form full FRF matrix, H [ns x no*ni]
%  H = [H11 H12      -> H = [H11 H21 H12 H22];
%       H21 H22 ];
H11 = A1(1,1)./(1j*w - root(1)) + A2(1,1)./(1j*w - root(2));
H21 = A1(2,1)./(1j*w - root(1)) + A2(2,1)./(1j*w - root(2));
H12 = A1(1,2)./(1j*w - root(1)) + A2(1,2)./(1j*w - root(2));
H22 = A1(2,2)./(1j*w - root(1)) + A2(2,2)./(1j*w - root(2));

% concat driving point info
H = [H11.' H21.' H12.' H22.']; % .' to preserve complex conj


%% FRF Plots using partial fraction expansion

in = 1;     % input dof (column index)
out = 1;    % output dof (row index)
tsize = 36;

% form FRF index, hInd
hInd = 1:no*ni;
hInd = reshape(hInd,no,ni);

% get individual contributions
h1 = A1(out,in)./(1j*w-root(1)); % mode 1 contribution
h2 = A2(out,in)./(1j*w-root(2)); % mode 2 contribution 


%% plots

figure; ah = axes;
ph1 = plot(w,mag2db(abs(H(:,hInd(out,in)))),'o-');
hold(ah,'all');
ph2 = plot(w,mag2db(abs(h1.')),'o-');
ph3 = plot(w,mag2db(abs(h2.')),'o-');

% mode 1 properties
set(ph1,'linewidth',2,...
            'markersize',2);
%         'markeredgecolor','b',...
%         'markerfacecolor','b',...

    
% mode 2 properties
set(ph2,'linewidth',2,...
        'markersize',2);
% mode 3 properties
set(ph3,'linewidth',2,...
        'markersize',2);
    
legend('Measurement', 'Mode: 1', 'Mode: 2');
grid(ah,'on'); grid(ah,'minor');
xlabel(ah,'Frequency [rad/sec]','fontsize',tsize,'fontname','Times New Roman');
ylabel(ah,'[db]' ,'fontsize',tsize,'fontname','Times New Roman');
set(ah,'fontsize',tsize,'fontname','Times New Roman');

% 
% % plot
% figure; ah = axes;
% plot(w,mag2db(abs([H(:,hInd(out,in)) h1.' h2.'])),'.-')




%% SVD for shapes

% loop freq for CMIF
uu = []; ss = []; vv = [];
for ii = 1:ns
    [ut,st,vt] = svd(reshape(H(ii,:),no,ni));
    uu(:,:,ii) = ut;
    ss(:,ii) = diag(st);
    vv(:,:,ii) = vt;
end

% plot
figure; ah = axes;
plot(w,mag2db(ss),'.-')
% legend({['H' num2str(out) num2str(in)]; 'Mode: 1'; 'Mode: 2';});
grid(ah,'on');
xlabel(ah,'Frequency [rad/sec]','fontsize',16,'fontname','Times New Roman');
ylabel(ah,'Singular Value' ,'fontsize',16,'fontname','Times New Roman');
set(ah,'fontsize',16,'fontname','Times New Roman');

% index poles
[~,peakLoc] = searchVector(w,[11.61 12.01]);

%% vector track

% artificial vector track
ff = 62;                            % freq index of mode switch
% ss = ss.';
s1 = [ss(1,1:ff-1) ss(2,ff:end)];   % mode 1
s2 = [ss(2,1:ff-1) ss(1,ff:end)];   % mode 2
SS = [s1; s2];


% plot
figure; ah = axes;
plot(w,mag2db(SS),'.-')
legend({'Mode: 1'; 'Mode: 2';});
grid(ah,'on');
xlabel(ah,'Frequency [rad/sec]','fontsize',16,'fontname','Times New Roman');
ylabel(ah,'Singular Value' ,'fontsize',16,'fontname','Times New Roman');
set(ah,'fontsize',16,'fontname','Times New Roman');


%% normalize mode shapes
% loop left singular vectors for mode shape
U = []; L = [];
for ii = 1:2    
    % normalize and save only real part
    U(:,ii) = real(normalizeMode(uu(:,1,peakLoc(ii))));
    
    % do not normalize, but save real part
    L(:,ii) = real(vv(:,1,peakLoc(ii)));
end

inLoc = 1:2; % index for input locations

% loop to scale left and right singular vectors
msf = []; Lm = [];
for ii = 1:2
    msf(ii) = pinv(L(:,ii))*U(inLoc,ii);
    Lm(:,ii) = msf(ii)*L(:,ii);
end


%% Curve Fitting - Modal Parameter Identification
%   use mode shapes to decouple MDOF to SDOF 
%   note: can be called DOF condensation
%         method used for eFRF

% loop spectral lines for eFRF
eH = [];
for ii = 1:ns
    % apply modal filter
    tt = U.'*reshape(H(ii,:),no,ni)*Lm;
    % save diagonal terms
    eH(ii,:) = diag(tt);
end


% pre-al
eHs = zeros(ns,ne); pols = zeros(ne,1); 
qr = zeros(ne,1);   ma = zeros(ne,1);
for ii = 1:2 % loop modes    
    % fit EFRF
    [ehs,pol,freqdamp,q,m] = eFRFfit2(eH.',U,Lm,w/2/pi,ii,peakLoc(ii),1:ni,5,4);    
    % save parameters
    eHs(:,ii) = ehs;    
    pols(ii) = pol;    
    frdmp(ii) = freqdamp;
    qr(ii) = q;         
    ma(ii) = m;
end

% plot eH and compare to true solution
figure; ah = axes;
plot(ah,w,mag2db(abs(eH)),'.-');
hold all
plot(ah,w,mag2db(abs(eHs)),'.-');
grid on


%% Modal Flexibility 

% dummy shape and scaling vars
UU = U;     % mode shape
Ma = ma;    % modal mass

% scale mode shapes w/ mass, and evaluate at jw = 0
flex = zeros(no,no);
for ii = 1:ne
    tmp = (UU(:,ii)*UU(:,ii).')/(-1*Ma(ii)*pols(ii));
    flex = flex + tmp + conj(tmp); % add each mode and complex conj
end


%% System Flexibility

%
Flex = zeros(no,no);
MA = [1/Q1 1/Q2];
for ii = 1:ne
    tmp = (V(:,ii)*V(:,ii).')/(-1*(MA(ii)*root(ii)));
    Flex = Flex + tmp + conj(tmp);
end

ff = [5; 0];    % force vector
Flex = inv(K);  % true flexibility matrix

DD = F*ff      % true disp vector
dd = flex*ff   % modal disp vector










%% Curve Fitting - True Solution (for comparison)
EH = [];
for ii = 1:ns
    % apply modal filter
    tt = V.'*reshape(H(ii,:),no,ni)*V;
    % save diagonal terms
    EH(ii,:) = diag(tt);
end

plot(ah,w,mag2db(abs(EH)),'.-');



%% plot shapes

xx = [0 l]; % x dimension
zz = U;    % modal amplitude

% plot
figure; ah = axes;
ph1 = plot(xx,zz(:,1),'o-');
hold(ah,'all');
ph2 = plot(xx,zz(:,2),'o-');
ph3 = plot(xx,zeros(size(xx)),'o-');
for kk = 1
% mode 1 properties
set(ph1,'linewidth',2,...
        'markeredgecolor','b',...
        'markerfacecolor','b',...
        'markersize',8);    
% mode 2 properties
set(ph2,'linewidth',2,...
        'color','r',...
        'markeredgecolor','r',...
        'markerfacecolor','r',...
        'markersize',8,...
        'displayname','Mode 2');
% un-deformed beam props
set(ph3,'linewidth',2,...
        'color','k',...
        'markeredgecolor','k',...
        'markerfacecolor','k',...
        'markersize',6);
% axes proprs
xlabel(ah,'Beam Length [in]','fontsize',18,'fontname','Times New Roman');
ylabel(ah,'Modal Amplitude' ,'fontsize',18,'fontname','Times New Roman');
set(ah,'fontsize',18,'fontname','Times New Roman');
ylim(ah,[-1 1]);
grid(ah,'on')
legend(ah,{'Mode 1' 'Mode 2'});
end













