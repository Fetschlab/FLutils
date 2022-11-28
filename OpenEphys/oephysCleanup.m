%% OpenEphys data cleanup
% after downloading from server, erase all the garbage we don't use.
% this will save disk space, subsequent loading time, and also make
% exploring the data easier.

% 10-19-20 CF
% must be run by getDataFromServer

% For now only erases the AUX data and chans 1-8 (which don't have data
% when using the 24-chan probe)


files = dir([localDir remoteFiles{n}]);
for f=3:length(files) % skip 1+2, they are are "." and ".."
    dot = strfind(files(f).name,'.');
    ext = files(f).name(dot(end)+1:end); % dot(end) because spikes files have two dots
    filename = fullfile(localDir,remoteFiles{n},files(f).name);
    
    if strcmp(ext,'continuous')
    
        % delete AUX channels
        if any(strfind(files(f).name,'AUX')) 
            disp(['deleting AUX chans: ' remoteFiles{n} '/' files(f).name])
            delete(filename);
        end
        
        % delete channels 1-8, not used currently
        ch = strfind(files(f).name,'CH');
        und = strfind(files(f).name,'_');
        if isempty(und)
            chan = str2double(files(f).name(ch+2:dot-1));
        else
            chan = str2double(files(f).name(ch+2:und(end)-1)); % und(end) because continuous files can have two underscores
        end
        if chan<9 
            disp(['deleting chans 1-8: ' remoteFiles{n} '/' files(f).name])
            delete(filename);
        end
        
    elseif strcmp(ext,'spikes')
        % delete spike files from electrodes 0-7, for the same reason
        zn = strfind(files(f).name,'0n');
        und = strfind(files(f).name,'_');
        if isempty(und)
            chan = str2double(files(f).name(zn+2:dot(end)-1)) + 1; % +1 because electrodes start with 0
        else
            chan = str2double(files(f).name(zn+2:und-1)) + 1; % +1 because electrodes start with 0            
        end
        if chan<9 
            disp(['deleting chans 1-8: ' remoteFiles{n} '/' files(f).name])
            delete(filename);
        end
        
    end
    
    
%     try % need try/catch to skip files like .DS_store
%         [~, ~, info] = load_open_ephys_data_faster(filename);
%         if strcmp(ext,'spikes')
%             
%         elseif strcmp(ext,'continuous')
%             chan = str2double(info.header.channel(3:end));
%         end
%     catch
%     end


end


