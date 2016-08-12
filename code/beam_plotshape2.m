function fh = beam_plotshape2(ah,y1,y2,out,L,nn)
%%
% 
% jdv 08202015; 11122015; 11142015

if isempty(ah)
    % create figure
    fh = figure; ah = axes; 
    fh.PaperPositionMode = 'auto';
    fh.Position = [50 250 1300 600];
end

% clear axes
cla(ah);

% gather spatial coordinates
xx = linspace(0,L/12,nn+2);       % x-dim in feet
yy = normalizeMode([y1(out) y2]); % normalize modes
yn = normalizeMode(y1);
xxo = xx(2:end-1);
xxo = xxo-xxo(1);
xxb = xxo(out);

% check sign
if sign(yy(1,1)) ~= sign(yy(1,2)) &&...
        sign(yy(2,1)) ~= sign(yy(2,2))
    % check second for roundoff
    if sign(yy(2,1)) ~= sign(yy(2,2))
        yy(:,2) = yy(:,2) * -1;
    end
end
if sign(yy(1)) ~= sign(yn(out(1))) &&...
        sign(yy(2)) ~= sign(yn(out(2)))
    yn = yn*-1;
end

% hold axes
hold(ah,'on');

% plot full resolution eigenvector
plot(ah,xxo,yn(:,1),'o-k',...
    'linewidth',2,...
    'markeredgecolor','none',...
    'markerfacecolor','k',...
    'markersize',14,...
    'displayname','Eigenvector');

% plot 1:1 eigenvector
plot(ah,xxb,yy(:,1),'o--b',...
    'linewidth',2,...
    'markeredgecolor','none',...
    'markerfacecolor','b',...
    'markersize',18,...
    'displayname','Exact');

% estimated
plot(ah,xxb,yy(:,2),'v--r',...
    'linewidth',3,...
    'markerfacecolor','r',...
    'markeredgecolor','none',...
    'markersize',18,...
    'displayname','Estimated');

% add legend before undeformed plot
legend(ah,'location','northeast');

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
