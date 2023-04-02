function h = plotmeanse(ymat,dim,xsamp,col,style,type)
%  plotmeanse(ymat,dim,xsamp,col,style)
%
% plots mean and error of ymat across dimension dim (default first non-singleton)
%
% xsamp - x-axis samples, should match plotting dimension of ymat
% col   - color (string, 3-length vector)
% String selection for how to plot SEM
%   'dash','shade','line','bar','none'

% Linewidths etc are hard coded for simplicity

% Adapted from Cavanagh et al 2016 eLife 
%
% Steven Jerjian March 2018

if ~ismatrix(ymat), error('ymat must be 2 dimensional matrix'); end


if nargin<2 || isempty(dim), dim = find(size(y)~=1,1); end
if nargin<4 || isempty(col), col = 'b'; end
if nargin<5 || isempty(style), style = 'shade'; end
if nargin<6 || isempty(type), type = 'sem'; end

validatestring(style,{'dash','shade','line','line0','bar','none'},'plotmeanse','style',5); %Jan2019

mn = nanmean(ymat,dim);
se = nanstd(ymat,[],dim)./sqrt(size(ymat,dim));

switch type
    case 'std'
        se = nanstd(ymat,[],dim);
    case 'prc'
        % not implemented
end

% transpose to row vector for consistency
if ~isrow(mn)
    mn=mn'; 
    se=se'; 
end

tf = ishold;
hold on;

if ischar(col),col=char2rgb(col); end

switch lower(style)
    
    case 'dash'
        h(1) = plot(xsamp,mn, 'LineWidth',1.5,'Color',col);
        h(2) = plot(xsamp,mn+se,'LineWidth',1,'Color',col,'Linestyle','--');
        h(3) = plot(xsamp,mn-se,'LineWidth',1,'Color',col,'Linestyle','--');
        
    case 'dot'
        h(1) = plot(xsamp,mn, 'LineWidth',1.5,'Color',col);
        h(2) = plot(xsamp,mn+se,'LineWidth',1,'Color',col,'Linestyle',':');
        h(3) = plot(xsamp,mn-se,'LineWidth',1,'Color',col,'Linestyle',':');
        
        
    case 'shade'
        h(1) = plot(xsamp,mn,'LineWidth',1.5,'Color',col);
        if ~any(isnan(mn)) 
            h(2) = patch([xsamp fliplr(xsamp)],[mn+se fliplr(mn-se)],col,'FaceAlpha',0.3,'EdgeAlpha',0);
        else % patch doesn't work if nans are in the data, this is a kluge
            naninds = find(isnan(mn));
            begs = [0 naninds(1:end)]+1;
            ends = naninds-1;
            if begs(end)<numel(mn),
                ends(end) = numel(mn);
            else
                begs = begs(begs<numel(mn));
            end
            for ii=1:length(begs)
                indsi = begs(ii):ends(ii);
                h(2) = patch([xsamp(indsi) fliplr(xsamp(indsi))],[mn(indsi)+se(indsi) fliplr(mn(indsi)-se(indsi))],col,'FaceAlpha',0.3,'EdgeAlpha',0);
            end
        end
        
    case 'line'
        h = errorbar(xsamp,mn,se,'LineWidth',1.5,'Color',col);
        
    case 'line0',
        h = errorbar(xsamp,mn,se,'LineWidth',1.5,'Color',col,'LineStyle','none');
        
    case 'bar'
        if size(col,1)==1
            [h,herr] = barwitherr(se,xsamp,mn);
%             set(h,'FaceColor',col,'LineWidth',1.5);
        
            % remove borders and set unidirectional error bars
%             set(h,'EdgeAlpha',0); set(herr, 'color',col,'LineWidth',1.5);
%             herr.YNegativeDelta(h.YData>=0) = NaN;
%             herr.YPositiveDelta(h.YData<0)  = NaN;
            
            % make borders and error bars black
            set(h,'EdgeColor','k'); set(herr, 'color','k','LineWidth',1.5);
            
            % 
            set(h,'FaceColor',col,'FaceAlpha',0.6,'LineWidth',1.5);
%             set(h,'EdgeColor',col);
%             set(herr, 'color',col,'LineWidth',1.5);
        else
            for c=1:numel(mn)
                h(c) = bar(xsamp(c),mn(c));
                set(h(c),'FaceColor',col(c,:),'LineWidth',1.5);
                set(h,'EdgeAlpha',0);
            end
        end
        
    case 'none'
        h = plot(xsamp,mn, 'LineWidth',1.5,'Color',col);
        
end

if ~tf,
    hold off;
end