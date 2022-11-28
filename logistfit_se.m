function [beta,llik,pred,se] = logistfit_se(data, guess)
% [beta,llik,pred,se] = logistfit_se(data)
% this function is similar to logist_fit.m but is slightly modified to
% measure the standard error of the bias directly.
% the old version of logistfit fits  data(:,end) = 1/(1+exp(-B*data(1:end-1))) which is
% equivalent to p = 1/(1+exp(-(B1+B2*C))). in the old version the bias, in units of
% coherence, is B1/B2. the old version does not provide the standard error of this ratio.
% logistfit_se, however, calculates this standard error. it fits p = 1/(1+exp(-B2*(B1+C))) 
% and returns the standard error for both B1 and B2
% more generally the fit is
%  data(:,end) = 1/(1+exp( -B(end)*(data(:,end-1)+B(1:end-1)*data(:,1:end-2)) ))
%
%  **** the function assumes that the column before the last in data contains the
%  coherence or its equivalent in the analysis ****
% 
% LOGISTFIT_SE  fits a logistic function to data using maximum likelihood
%	maximization under binomial assumptions. Data must be organized as rows
%	of observations. The 1st column should be a 1 if you want to estimate a
%	'bias' term, columns 2 thru n-1 are the values for the independent
%	variables, with column n-1 containing the coherence or its equivalent.
%	The last column contains 1 or 0 for the choice. 
%
%	usage: [beta,llik,pred,se] = logistfit(data)
%	Return values are
%	   BETA, the coefficients corresponding to columns 1:m-1 of data.
%	   LLIK, the log likelihood of the fit
%	   PRED, a vector representing the prob(y=1) under the fit.
%     SE,   a vector of the standard errors of the estimated parameters
%
% created by RK, 05/16/2007

n = size(data,2);

% generate guess
if nargin < 2 || isempty(guess),
    % to generate a 1st guess, we will simply convert all the 1's to .99 and
    % all the 0's to .01's and do a quick regression on the logits.  This is
    % roughly converting the 1's to logit of 3 and 0's to logit of -3
    y = (6 * data(:,n)) - 3;
    q = data(:,1:n-1)\y;			
    guess = [q(1:end-1)/q(end); q(end)];       % this is a 1st guess
end;

if (nargout>3)
    % standard errors
    % The covariance matrix is the negative of the inverse of the hessian of the natural logarithm
    % of the probability of observing the data set given the optimal parameters.
    options = optimset ( 'MaxFunEvals', 400*length(guess), 'MaxIter', 400*length(guess), ...
                         'LargeScale' , 'off' , 'Display' , 'off' , 'GradObj' , 'off' , 'Hessian' , 'off' );
    [beta,fval,~,~,~,hessian] = fminunc(@(x) logiterr_se(x,data),guess,options);
    se = sqrt(diag(inv(hessian)));
else
    options = optimset('MaxFunEvals',400*length(guess),'MaxIter',400*length(guess),'Display','off');
    [beta,fval] = fminsearch(@(x) logiterr_se(x,data),guess,options);      
end;

if nargout > 1
  llik = -fval;
end

if nargout > 2
  % calculate the predicted prob(resp=1) values for each observation.
  qcol = beta(:);			% we need a column vector of params
  z = exp( (data(:,1:n-2)*qcol(1:end-1)+data(:,n-1))*qcol(end) );
  pred = z ./ (1 + z); 			% calculate the prob of getting a 1
end


% data = ones, I_s, coh*Is, coh, choice


function err = logiterr_se(q,data)

[~,n] = size(data);
qcol = q(:);				% we need a column vector of params
z = exp( (data(:,1:n-2)*qcol(1:end-1)+data(:,n-1))*qcol(end) );

%   that is: [all data columns except coh, multiplied by their respective beta terms, PLUS
%             the coherence itself] TIMES the coherence beta term
%   aka,:
%   B_c*[B1 + B2*Is + B3*Is*C + C]

p1 = z ./ (1 + z);			% calculate the prob of getting a 1
% for each data entry.

% Now calculate the joint probability of obtaining the data set conceived as
% a list of Bernoulli trials.  This is just p1 for trials = 1 and 1-p1 for
% trials of 0.

p = [p1(data(:,n)==1); 1-p1(data(:,n)==0)];
err = -sum(log(p));			% minus log likelihood




