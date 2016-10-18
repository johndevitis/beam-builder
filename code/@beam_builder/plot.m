function fh = plot(self,fh)
%% plot
% 
% 
% 
% author: john devitis
% create date: 18-Oct-2016 16:13:20
	
	if nargin < 2
        % create figure
        fh = figure; 
        ah = axes; 
        fh.PaperPositionMode = 'auto';
        fh.Position = [100 100 1300 600];
    else
        ah = get(fh,'Children');
    end
	hold(ah,'on');
    
    zs = zeros(size(self.x)); % zeros for undeformed
    % plot all nodes - undeformed beam
	plot(ah,self.x,zs,'o:k',...
        'linewidth',2,...
        'markerfacecolor','none',...
        'markersize',11);
    
    hold(ah,'off');
end
