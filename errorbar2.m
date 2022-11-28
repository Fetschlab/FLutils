function [hsym,hxe,hye] = errorbar2(x,y,xe,ye,markertype,linewidth)
% plots a double errorbar plot.
% [hsym,hxe,hye] = errorbar2(x,y,xe,ye)
% Routine plots x+-xe vs. y+-ye
%
% input args
% ~~~~~~~~~~
% x  vector of x values
% y  vector of y values
% xe  vector of x error 
% ye  vector of y errors
% 
% return args
% ~~~~~~~~~~~
% hsym   handle to symbols
% hex    handle to plot of the x error lines
% hey    handle to plot of the y error lines

% 3/17/98 mns wrote it
% 7/28/03 replaced nans with repmat(nan,...)
% 3/18/04  hold state now respected

if nargin<5
    markertype = 'o';
end

if nargin<6
    linewidth = 1; % CRF
end

x = x(:);
y = y(:);
xe = xe(:);
ye = ye(:);

hold_state = ishold;  %test the status of hold

% This is a trick to get line segments
qy = [y y nan(size(y))]';
qx = [x x nan(size(x))]';
qye = [y-ye y+ye nan(size(y))]';
qxe = [x-xe x+xe nan(size(x))]';

hxe = plot(qxe(:),qy(:),'k-'); hold on;
hye = plot(qx(:),qye(:),'k-'); hold on;

% Plot the symbols last, so they are in front.
hsym = plot(x,y,markertype,'LineWidth',linewidth);
% set(hsym,'MarkerFaceColor','w');

if hold_state
    hold on;
else
    hold off;
end



