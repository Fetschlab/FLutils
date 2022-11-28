function f=printfig(f,name,filetype,siz,append)
% printfig for publication
% requires export_fig package for some uses
%
% f: handle to figure, defaults to current figure
% name: name to save figure as, as string, defaults to newfig
% filetype: ps, eps, png, pdf, defaults to pdf
% siz: desired [x,y] size of figure in centimeters, defaults to current
% append: flag whether to append figure to already saved file (if file
% doesn't exist it will be created)

if nargin < 5, append = 0; end
if nargin < 3, filetype = 'pdf'; end
if nargin < 2, name = 'newfig'; end

if nargin < 1 ||isempty(f), f = gcf; end

set(f,'Units','centimeters');
if nargin < 4 ||isempty(siz), siz=f.Position(3:4); end

x=siz(1);
y=siz(2);

% we want figure on paper to look exactly as it does in Matlab
set(f,'PaperUnits','centimeters');
set(f,'PaperPositionMode', 'manual');
set(f,'PaperSize',[x,y],'PaperPosition',[0,0,x,y])
set(f,'renderer', 'painters');

%to see on screen as on print
set(f,'position',[10,10,x,y])



if append
    switch filetype
        case 'ps' , print(f,'-dpsc2',name,'-append');
        case 'eps', print(f,'-depsc2',name,'-append');
        case 'png', export_fig(name,f,'-png','-nocrop','-append');
        case 'pdf', export_fig(name,f,'-pdf','-nocrop','-append');
    end
else
    switch filetype
        case 'ps' , print(f,'-dpsc2',name); 
        case 'eps', print(f,'-depsc2',name);
        case 'png', export_fig(name,f,'-png','-nocrop');
        case 'pdf', print(f,'-dpdf',name);
    end 
end