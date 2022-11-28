% grabs OpenEphys data from NAS as specified, then loads desired variables
% from .PDS files into a simpler data structure

% 10-18-20

% requires password(s) to be entered if a auth key not available 
% see e.g.:
% https://www.howtogeek.com/66776/how-to-remotely-copy-files-over-ssh-without-entering-your-password/

clear all

%% decide which files to load

subject = 'hanzo';
paradigm = 'Dots';

dateRange = 20201102:20201102;


% dateRange = 20201001:20201031;
% dateRange = 20201101:20201130;
% dateRange = 20201201:20201231;
% dateRange = 20210101:20210131;
% dateRange = 20210201:20210228;
% dateRange = 20210301:20210331;


% dateRange = 20210101:20210309;


dateStr = num2str(dateRange(1));
for d = 2:length(dateRange)
    dateStr = [dateStr newline num2str(dateRange(d))];
end

localDir = ['/Users/chris/Documents/MATLAB/OpenEphys_data/' subject '/'];
remoteDir = ['/var/services/homes/fetschlab/data/' subject '_neuro/'];


%% get PDS files from server
% will skip files that already exist locally, unless overwrite set to 1

overwriteLocalFiles = 0; % set to 1 to always use the server copy

tic
getDataFromServer_sortedOnly % searches for .psort files and only downloads them and their associated chans/runs
toc

% save currentFileList so you don't have to repeat the time consuming step
if exist('currentFilename','var')

if sum(diff(dateRange)>1)==0
    file = [subject '_neuro_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];
elseif sum(diff(dateRange)>1)==1
    file = [subject '_neuro_' num2str(dateRange(1)) '-' num2str(dateRange(diff(dateRange)>1)) '+' num2str(dateRange(find(diff(dateRange)>1)+1)) '-' num2str(dateRange(end)) '.mat'];
else
    warning('Multiple breaks in date range, could overwrite existing file with different intervening dates!');
    disp('Press a key to continue');
    pause
    file = [subject '_neuro_' num2str(dateRange(1)) '---' num2str(dateRange(end)) '.mat'];
end

save([localDir(1:length(localDir)-length(subject)-1) file], 'data');

end




%% data browser (work in progress, may be abandoned)

% % oephysDataBrowser






