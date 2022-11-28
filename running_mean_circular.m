function [x_,y_,se_] = running_mean_circular(x, y, chunklen, flag, data_type)
%
% function [x_,y_, se_] = running_mean ( x , y , chunklen , flag )
%
% calculate the running mean of y on the axis defined by x
%
% if flag is 'data_point', the default, chunklen defines the
% number of data points, around each value of x, that will be used for
% the calculation of running mean.
% if flag is 'window_width' chunklen will be interpreted as the size of the
% window, in units of x, that will be used for the calculation of running 
% mean.
% 

dbstop if error

% RK, 12/19/2006 
% RK, 2/24/2010, the function also returns se_, the standard error of the running mean
% CF, 9/23/2021, modified for circular data; x must be angle in degrees

if nargin<4 || isempty(flag)
    flag = 'data_point';
end

if nargin<5 || isempty(data_type)
    data_type = 'continuous';
end

I = ~isnan(x) & ~isnan(y);
x = x(I);
y = y(I);

[x_,I] = sort(x);
y = y(I);
y_ = nan(size(x_));
se_ = nan(size(y_));

if chunklen > length(y)
    error('chunklen is longer than y; cannot wrap twice!')
end

switch flag
    case 'data_point'        
        for d = 1 : length(x_)
%             I = find(x_==x_(d)); % CF: don't understand this. why does it matter if there are duplicate x vals?
%             prop = round((chunklen-length(I))/2);
%             if prop<0
%                 prop = 0;
%             end
%             minI = I(1)-prop;
%             maxI = I(end)+prop;

            minI = d - floor(chunklen/2);
            maxI = d + floor(chunklen/2);

            if minI<0 && maxI>length(y)
                error('window appears to be longer than y; cannot wrap twice!')
            end
            
            if minI<=0
                set1 = y(end+minI:end);
                set2 = y(1:maxI);
            elseif maxI>length(y)
                set1 = y(minI:end);
                set2 = y(1:maxI-length(y));
            else
                set1 = y(minI:maxI);
                set2 = nan;
            end
            
            try y_(d) = nanmean([set1 set2]); catch; y_(d) = nanmean([set1;set2]); end

            if nargout>2
                switch data_type
                    case 'continuous'
                        try se_(d) = nanstd([set1 set2])/sqrt(maxI-minI+1); catch; se_(d) = nanstd([set1;set2])/sqrt(maxI-minI+1); end
                    case 'binary'
                        se_(d) = sqrt(y_(d)*(1-y_(d))/(maxI-minI+1));
                end
            end
        end
        
    case 'window_width'
        for d = 1 : length(x_)
            minV = x_(d)-chunklen/2;
            maxV = x_(d)+chunklen/2;
            if minV < 0 && maxV > x_(end)
                error('window appears to be longer than y; cannot wrap twice!')
            end
            if minV < 0
                I = x_>360+minV | x_<maxV;
            elseif maxV > 360
                I = x_>minV | x_<maxV-360;
            else
                I = x_>minV & x_<maxV;
            end
            y_(d) = mean(y(I));
            if nargout>2
                switch data_type
                    case 'continuous'
                        se_(d) = std(y(I))/sqrt(sum(I));
                    case 'binary'
                        se_(d) = sqrt(y_(d)*(1-y_(d))/sum(I));
                end
            end
        end        
end

    

    
