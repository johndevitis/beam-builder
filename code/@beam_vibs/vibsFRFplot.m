function fh = vibsFRFplot(Hs,hh,in,out,w)
%% Vibs Example - Construct FRF - Mode Contributions
%
% logmag plot of Hpq with dashed mode contributions
% complex conjugate solution
%
% jdv 06232015; 07232015; 08162015; 11112015

% plot frf
fprintf('H%d%d Magnitude\n',out,in);
    
% create figure
fh = figure; 
ah = axes; 
fh.PaperPositionMode = 'auto';
fh.Position = [100 100 1300 600];

% plot frf
plot(w,mag2db(abs(Hs)),'o-k','linewidth',2,'markersize',2);
hold all
plot(ah,w,mag2db(abs(hh)),'--','linewidth',1.5);
hold off

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
% ybnds = [-110 -50]; 
xres = 10;
yres = 10;
xlabel(ah,'Frequency [rad/sec]');
set(ah,'fontsize',18,'fontname','times new roman');
% set(ah,'xtick',xbnds(1):xres:xbnds(2));
% set(ah,'ytick',ybnds(1):yres:ybnds(2));
% xlim(ah,xbnds); 
ylabel(ah,'[db]');
% ylim(ah,ybnds);



