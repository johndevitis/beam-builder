function plotShape(self,bm,mode)
%% plotShape
% 
% bm.
%   x -  longitudinal coordinates
%   bc - boundary condition indices
% 
% author: john devitis
% create date: 18-Oct-2016 18:12:13
	
fh = bm.plot();
ah = get(fh,'Children');
hold(ah,'on');
   
xx = bm.x;
yy = zeros(size(xx));

xxo = bm.x(2:end-1);
yyo = self.V(:,mode);
    
% plot full resolution eigenvector
plot(ah,xxo,yyo,'o-k',...
    'linewidth',2,...
    'markeredgecolor','none',...
    'markerfacecolor','k',...
    'markersize',14,...
    'displayname','Eigenvector');
	
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
end
