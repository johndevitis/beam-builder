function ah = beam_plotshape(ah,V,mode,L,nn)
%%
% 
% jdv 08202015; 11142015

if isempty(ah)
    % create figure
    fh = figure; 
    ah = axes; 
    fh.PaperPositionMode = 'auto';
    fh.Position = [100 100 1300 600];
end

% gather spatial coordinates
xx = linspace(0,L/12,nn+2);    % x-dim in feet
yy = normalizeMode(V);         % normalize modes
xxo = xx(2:end-1);
xxo = xxo-xxo(1);

% plot full resolution eigenvector
plot(ah,xxo,yy(:,1),'o-k',...
    'linewidth',2,...
    'markeredgecolor','none',...
    'markerfacecolor','k',...
    'markersize',14,...
    'displayname','Eigenvector');

hold(ah,'all');

% check for second mode
if length(mode) == 2
    % estimated
    plot(ah,xxo,yy(:,2),'v--r',...
        'linewidth',3,...
        'markerfacecolor','r',...
        'markeredgecolor','none',...
        'markersize',18,...
        'displayname','Estimated');
end
    
% plot un-deformed
plot(xxo,zeros(size(xxo)),'o:k',...
    'linewidth',2,...
    'markerfacecolor','none',...
    'markersize',11);

hold(ah,'off');

% format
tsize = 26; % text size
xlabel(ah,'Beam Length [ft]');
ylabel(ah,'Modal Amplitude');
set(ah,'fontsize',tsize,'fontname','Times New Roman');
ylim(ah,[-1.1 1.1]);
xlim(ah,[min(xxo)-5 max(xxo)+5])
grid(ah,'on')
grid(ah,'minor')
box on
