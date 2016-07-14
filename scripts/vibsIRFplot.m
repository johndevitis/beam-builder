function fh = vibsIRFplot(hs,h,in,out,fs,l)
%% Vibs Example - Impulse Resopnse Function
%
%
% jdv 06232015; 07232015; 08162015; 11112015

% time sampling parameters
dt = 1/fs;
t  = 0:dt:l-dt;

% create figure
fh = figure; 
ah = axes; 
fh.PaperPositionMode = 'auto';
fh.Position = [200 200 1300 600];

plot(t,hs,'o-k','linewidth',2,'markersize',2);
hold(ah,'all');
plot(t,h,'--','linewidth',1.5);
hold(ah,'off');

% format plot
grid(ah,'on'); grid(ah,'minor');
xlabel(ah,'Time [sec]');
ylabel(ah,'Displacement [in]')
set(ah,'fontname','times new roman');
set(ah,'fontsize',18);

% form legend
ne = size(h,2);
lg = {['H' num2str(out) num2str(in)]};
for ii = 2:ne+1
    lg{ii} = ['Mode: ' num2str(ii-1)];
end
lh = legend(lg,'location','northeast');

