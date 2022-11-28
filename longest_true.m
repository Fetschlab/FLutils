function [longest,st,lens] = longest_true(n)
% [longest,st,lens] = longest_true(n)
%
%
% find longest sequence of ones in a vector 
%
% OUTPUTS:
% longest - length of longest sequence of consecutive true values
% st      - indices of starts of each sequence of true values
% lens    - length of all clusters

%%
% deal with empty or all zeros case

if ~isrow(n)
    n = n';
end

if isempty(nonzeros(n))
    longest = 0; st=[]; lens=[]; return
end

% deal with cell input
if iscell(n)
    longest = max(cellfun(@numel,regexp(strrep(n,' ',''),'1+','match'))); 
    st=[]; lens=[];
    return
end

if ischar(n)
    n = str2num(n); %#ok<ST2NM>
end

% if nargout>1
    n0      = (cumsum(~n)+1) .* n;
    idx     = find(n0);
    st      = idx([0 find(diff(idx)>1)]+1);            % starts of clusters
%     en      = idx([find(diff(idx)>1) numel(idx)]);   % ends of clusters
    lens    = nonzeros(accumarray(nonzeros(n0), 1));   % lengths of clusters
    longest = max(lens);
    
% else
%     longest = max( accumarray ( nonzeros((cumsum(~n)+1) .* n), 1 ));
%     st = [];
%     lens = [];
% end


