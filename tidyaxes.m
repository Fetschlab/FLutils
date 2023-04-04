function ax=tidyaxes(ax,fontsize)

if nargin < 1, ax = gca; end
if nargin < 2, fontsize = 14; end
    
if nargin==2
    set(ax,'FontSize',fontsize);
    set(get(ax,'XLabel'),'FontSize',fontsize);
    set(get(ax,'YLabel'),'FontSize',fontsize);
    set(get(ax,'Title'),'FontSize',fontsize+2);
end
set(get(ax,'Title'),'FontWeight','bold');

ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
ax.TickDir         = 'out';

ax.TickLength = [0.02 0.025];

ax.Color           = [1 1 1];
ax.Box             = 'off';
