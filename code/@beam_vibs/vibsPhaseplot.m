function fh = vibsPhaseplot(Hs,hh,in,out,w)
%% Vibs Example - Construct Phase - Mode Contributions
%
% phase plot of Hpq with dashed mode contributions
% complex conjugate solution
%
% jdv 06232015; 07232015; 08162015; 11112015

% create figure
fh = figure; 
ah = axes; 
fh.PaperPositionMode = 'auto';
fh.Position = [150 150 1300 600];

% plot phase
plot(ah,w,angle(Hs)*180/pi,'o-k','linewidth',2,'markersize',2);
hold(ah,'all');
plot(ah,w,angle(hh)*180/pi,'--','linewidth',1.5);

% form legend
ne = size(hh,2);
lg = {['H' num2str(out) num2str(in)]};
for ii = 2:ne+1
    lg{ii} = ['Mode: ' num2str(ii-1)];
end
legend(lg,'location','northeast');

% format
grid(ah,'on'); 
grid(ah,'minor');
% xbnds = [0 w(end)];
ybnds = [-190 190]; 
xres = 10;
yres = 10;
xlabel(ah,'Frequency [rad/sec]');
set(ah,'fontsize',18,'fontname','times new roman');
% set(ah,'xtick',xbnds(1):xres:xbnds(2));
set(ah,'ytick',ybnds(1):yres:ybnds(2));
% xlim(ah,xbnds); 
ylabel(ah,'Phase [deg]');
% ylim(ah,ybnds);
set(ah,'ytick',-180:45:180);
    
end