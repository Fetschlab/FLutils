function [data, timestamps, info] = load_open_ephys_data_faster(filename, varargin)

%
% [data, timestamps, info] = load_open_ephys_data(filename, [outputFormat])
%
%   Loads continuous, event, or spike data files into MATLAB.
%
%   Inputs:
%
%     filename: path to file
%     outputFormat: (optional) If omitted, continuous data is output in 
%                   double format and is scaled to reflect microvolts. If
%                   this argument is 'unscaledInt16' and the file contains
%                   continuous data, the output data will be in int16
%                   format and will not be scaled; this data must be
%                   manually converted to a floating-point format and
%                   multiplied by info.header.bitVolts to obtain microvolt
%                   values. This feature is intended to save memory for
%                   operations involving large amounts of data.
%
%
%   Outputs:
%
%     data: either an array continuous samples (in microvolts unless outputFormat is specified, see above),
%           a matrix of spike waveforms (in microvolts),
%           or an array of event channels (integers)
%           OR a struct with trial events (messages) sent over network -CF
%
%     timestamps: in seconds
%
%     info: structure with header and other information
%
%
%
%   DISCLAIMER:
%
%   Both the Open Ephys data format and this m-file are works in progress.
%   There's no guarantee that they will preserve the integrity of your
%   data. They will both be updated rather frequently, so try to use the
%   most recent version of this file, if possible.
%
%

%
%     ------------------------------------------------------------------
%
%     Copyright (C) 2014 Open Ephys
%
%     ------------------------------------------------------------------
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     <http://www.gnu.org/licenses/>.
%


[~,~,filetype] = fileparts(filename);
if ~any(strcmp(filetype,{'.events','.continuous','.spikes'}))
    error('File extension not recognized. Please use a ''.continuous'', ''.spikes'', or ''.events'' file.');
end

bInt16Out = false;
if nargin > 2
    error('Too many input arguments.');
elseif nargin == 2
    if strcmpi(varargin{1}, 'unscaledInt16')
        bInt16Out = true;
    else
        error('Unrecognized output format.');
    end
end

fid = fopen(filename);
fseek(fid,0,'eof');
filesize = ftell(fid);


% CF: for ordinary event files, the header contains evaluatable commands
% which generate the struct 'info'. This is not true for our network events
% files ("messages"), so we skip this in that case, and run a separate
% loading routine in the second switch-case stack below
if any(strfind(filename,'messages'))
    NUM_HEADER_BYTES = 0;
    info = [];
    version = 0;
else
    NUM_HEADER_BYTES = 1024;
    fseek(fid,0,'bof');
    hdr = fread(fid, NUM_HEADER_BYTES, 'char*1');
    info = getHeader(hdr);
    if isfield(info.header, 'version')
        version = info.header.version;
    else
        version = 0.0;
    end
end


switch filetype
    case '.events'
        bStr = {'timestamps' 'sampleNum' 'eventType' 'nodeId' 'eventId' 'data' 'recNum'};
        bTypes = {'int64' 'uint16' 'uint8' 'uint8' 'uint8' 'uint8' 'uint16'};      
        bRepeat = {1 1 1 1 1 1 1};
        dblock = struct('Repeat',bRepeat,'Types', bTypes,'Str',bStr);
        if version < 0.2, dblock(7) = [];  end
        if version < 0.1, dblock(1).Types = 'uint64'; end
    case '.continuous'
        SAMPLES_PER_RECORD = 1024;
        bStr = {'ts' 'nsamples' 'recNum' 'data' 'recordMarker'};
        bTypes = {'int64' 'uint16' 'uint16' 'int16' 'uint8'};
        bRepeat = {1 1 1 SAMPLES_PER_RECORD 10};
        dblock = struct('Repeat',bRepeat,'Types', bTypes,'Str',bStr);
        if version < 0.2, dblock(3) = []; end
        if version < 0.1, dblock(1).Types = 'uint64'; dblock(2).Types = 'int16'; end
    case '.spikes'
        num_channels = info.header.num_channels;
        num_samples = 40; 
        bStr = {'eventType' 'timestamps' 'timestamps_software' 'source' 'nChannels' 'nSamples' 'sortedId' 'electrodeID' 'channel' 'color' 'pcProj' 'samplingFrequencyHz' 'data' 'gain' 'threshold' 'recordingNumber'};
        bTypes = {'uint8' 'int64' 'int64' 'uint16' 'uint16' 'uint16' 'uint16' 'uint16' 'uint16' 'uint8' 'float32' 'uint16' 'uint16' 'float32' 'uint16' 'uint16'};
        bRepeat = {1 1 1 1 1 1 1 1 1 3 2 1 num_channels*num_samples num_channels num_channels 1};
        dblock = struct('Repeat',bRepeat,'Types', bTypes,'Str',bStr);
        if version < 0.4,  dblock(7:12) = []; dblock(8).Types = 'uint16'; end
        if version == 0.3, dblock = [dblock(1), struct('Repeat',1,'Types','uint32','Str','ts'), dblock(2:end)]; end
        if version < 0.3, dblock(2) = []; end
        if version < 0.2, dblock(9) = []; end
        if version < 0.1, dblock(2).Types = 'uint64'; end
end
blockBytes = str2double(regexp({dblock.Types},'\d{1,2}$','match', 'once')) ./8 .* cell2mat({dblock.Repeat});
numIdx = floor((filesize - NUM_HEADER_BYTES)/sum(blockBytes));

switch filetype
    case '.events'
        if any(strfind(filename,'messages'))
            % CF: read network events (messages)
                % utterly mystified as to why fid does not persist into my
                % new readMessages function (below), nor why fread fails
                % unless I do a new fopen.. but cest la vie.
            fid = fopen(filename);
            raw = fread(fid, inf, 'char*1');
            data = readMessages(raw);
            timestamps = [];
        else
            timestamps = segRead('timestamps')./info.header.sampleRate;
            info.sampleNum = segRead('sampleNum');
            info.eventType = segRead('eventType');
            info.nodeId = segRead('nodeId');
            info.eventId = segRead('eventId');
            data = segRead('data');
            if version >= 0.2, info.recNum = segRead('recNum'); end
        end
    case '.continuous'
        if nargout>1
            info.ts = segRead('ts');
        end            
        info.nsamples = segRead('nsamples');
        if ~all(info.nsamples == SAMPLES_PER_RECORD)&& version >= 0.1, error('Found corrupted record'); end
        if version >= 0.2, info.recNum = segRead('recNum'); end
        
        % read in continuous data
        if bInt16Out
            data = segRead_int16('data', 'b');
        else
            data = segRead('data', 'b') .* info.header.bitVolts;
        end
        
        if nargout>1 % do not create timestamp arrays unless they are requested
            timestamps = nan(size(data));
            current_sample = 0;
            for record = 1:length(info.ts)
                timestamps(current_sample+1:current_sample+info.nsamples(record)) = info.ts(record):info.ts(record)+info.nsamples(record)-1;
                current_sample = current_sample + info.nsamples(record);
            end
            timestamps = timestamps./info.header.sampleRate;
        end
    case '.spikes'
        timestamps = segRead('timestamps')./info.header.sampleRate;
        info.source = segRead('source');
        info.samplenum = segRead('nSamples');
        info.gain = permute(reshape(segRead('gain'), num_channels, numIdx), [2 1]);
        info.thresh = permute(reshape(segRead('threshold'), num_channels, numIdx), [2 1]);
        if version >= 0.4, info.sortedId = segRead('sortedId'); end
        if version >= 0.2, info.recNum = segRead('recordingNumber'); end
        data = permute(reshape(segRead('data'), num_samples, num_channels, numIdx), [3 1 2]);
        data = (data-32768)./ permute(repmat(info.gain/1000,[1 1 num_samples]), [1 3 2]);
end
fclose(fid);



function seg = segRead_int16(segName, mf)
    %% This function is specifically for reading continuous data. 
    %  It keeps the data in int16 precision, which can drastically decrease
    %  memory consumption
    if nargin == 1, mf = 'l'; end
    segNum = find(strcmp({dblock.Str},segName));
    fseek(fid, sum(blockBytes(1:segNum-1))+NUM_HEADER_BYTES, 'bof'); 
    seg = fread(fid, numIdx*dblock(segNum).Repeat, [sprintf('%d*%s', ...
        dblock(segNum).Repeat,dblock(segNum).Types) '=>int16'], sum(blockBytes) - blockBytes(segNum), mf);
    
end

function seg = segRead(segName, mf)
    if nargin == 1, mf = 'l'; end
    segNum = find(strcmp({dblock.Str},segName));
    fseek(fid, sum(blockBytes(1:segNum-1))+NUM_HEADER_BYTES, 'bof'); 
    seg = fread(fid, numIdx*dblock(segNum).Repeat, sprintf('%d*%s', ...
        dblock(segNum).Repeat,dblock(segNum).Types), sum(blockBytes) - blockBytes(segNum), mf);
    
end

end
function info = getHeader(hdr)
    eval(char(hdr'));
    info.header = header;
end


% CF, based on our early modified version of load_open_ephys_data
function events = readMessages(raw)
    
    str = char(raw');
    char10 = [1; find(raw==10)]; % 'carriage return' character
    edata = cell(length(char10)-1,1);
    T = 1; % trial counter
    D = 1; % counter for direction
    for n = 1:length(char10)-1
        edata{n} = [str(char10(n):char10(n+1)-1) ' ']; % add a space at the end so sp(2) always exists
        sp = strfind(edata{n},' ');
        thisStr = edata{n}(sp(1)+1:sp(2)-1);
        if strcmp(thisStr,'TrialStart') %mvl: Find all TrialStart cells
            trInd(T,1) = n; %mvl: index of all "TrialStart"
            T = T+1;
        end
        if strcmp(thisStr, 'Direction')
            trInd_Direction(D,1) = n; %mvl: index of all "Direction"
            D = D+1;
        end
    end
    
    % initialize
    events.TrialNum = nan(T-1,1);
    timedEvents = {'TrialStart','MotionStart','MotionEnd','Breakfix','GoodTrial'};
    trialParams = {'Coherence','Duration','GoodOrBadTrial','Speed','Diameter'};
    both = {'ApertureX','ApertureY','Direction'};
    for k = 1:length(timedEvents)
        eval(['events.' timedEvents{k} '=nan(T-1,1);']);
    end
    for k = 1:length(trialParams)
        eval(['events.' trialParams{k} '=nan(T-1,1);']);
    end
    for k = 1:length(both)
        eval(['events.' both{k} '=cell(T-1,1);']);
    end

    trInd(end+1) = length(char10); % add an extra dummy index for end of last trial
    for t = 1:T-1 % loop over trials
        % the first event in a trial is TrialStart, by definition, and this
        % line also contains the trial number, which may differ from the
        % index n, so we should keep track of it:
        sp = strfind(edata{trInd(t)}, ' '); %find spaces, which should seperate the Time of Event, Even Name, and Any other Info (ex. Coherence Value)
        events.TrialNum(t) = str2double(edata{trInd(t)}(sp(2)+1:end-1)); %mvl:get trial number

        for n = trInd(t):trInd(t+1)-1 % for each trial, loop over events (all events within each 'TrialStart')
            sp = strfind(edata{n},' ');
            thisStr = edata{n}(sp(1)+1:sp(2)-1);
            if ismember(thisStr,timedEvents) %mvl% Is the event message a TimedEvent or Parameter Event
                ts = str2double(edata{n}(1:sp(1)-1)); %mvl: Everything before the first space (which should contain the time)
                eval(['events.' thisStr '(t)=ts;']);
            elseif ismember(thisStr,trialParams) %mvl: whatever isnt fill is left off as NaN, which is convenient
                val = str2double(edata{n}(sp(2)+1:end-1));
                eval(['events.' thisStr '(t)=val;']);
            elseif ismember(thisStr,both)
                ts = str2double(edata{n}(1:sp(1)-1));
                val = str2double(edata{n}(sp(2)+1:end-1));
                eval(['events.' thisStr '{t}(end+1,1)=ts;']); % store the ts and val together for events that occur multiple times per trial
                eval(['events.' thisStr '{t}(end,2)=val;']); % (the previous end+1 takes care of inserting a new line, so this gets just 'end')
            end % otherwise skip it (Processor, Software, IsRecording, etc)
        end
        
    end

% % %     % MVL's separate loop for Direction, may be obsolete
% % %     trInd_Direction(end+1) = length(char10);
% % %     for t = 1:D-1 % loop over trials
% % %         % the first event in a trial is TrialStart, by definition, and this
% % %         % line also contains the trial number, which may differ from the
% % %         % index n, so we should keep track of it:
% % %         sp = strfind(edata{trInd_Direction(t)}, ' '); %find spaces, which should seperate the Time of Event, Even Name, and Any other Info (ex. Coherence Value)
% % %         %events.TrialNum(t) = str2double(edata{trInd_Direction(t)}(sp(2)+1:end-1)); %mvl:get trial number
% % %         for n = trInd_Direction(t):trInd_Direction(t+1)-1 % for each trial, loop over events (all events within each 'TrialStart')
% % %             sp = strfind(edata{n},' ');
% % %             thisStr = edata{n}(sp(1)+1:sp(2)-1);
% % %             if ismember(thisStr,'Direction') %mvl% Is the event message a TimedEvent or Parameter Event
% % %                 ts = str2double(edata{n}(1:sp(1)-1)); %mvl: Everything before the first space (which should contain the time)
% % %                 eval(['events.T' thisStr '(t)=ts;']);
% % %             end
% % %             if ismember(thisStr,'Direction') %mvl: whatever isnt fill is left off as NaN, which is convinient
% % %                 val = str2double(edata{n}(sp(2)+1:end-1));
% % %                 eval(['events.' thisStr '(t)=val;']);
% % %             end % otherwise skip it (Processor, Software, IsRecording, etc)
% % %         end
% % %     end


end
