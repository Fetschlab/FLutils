function x = naninterp(x,method)
% x = NANINTERP(x,method)
% find NaN values of underlying function (i.e. non-NaN values at non-NaN timepoints)
% wrapper for MATLAB's interp1

if nargin < 2,
    method = 'linear';
end

x(isnan(x)) = interp1(find(~isnan(x)),x(~isnan(x)),find(isnan(x)),method);