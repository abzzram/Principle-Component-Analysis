% PARAMETERIZE DATA
%   Takes single cell traces and breaks them into scalar parameters
%   outputs a nCell x Parameter matrix for use in PLSR
% -----------------------------------------------------------------------
% INPUTS
% d - single cell trace data
%   Accepted formats: cell array of traces (from ct_dataproc), cell array
%   containing multiple channel vaclubes, or cell array of
%   nCell x Time matrices
%   channel
%   *** IF YOUR DATA HAS MULTIPLE CHANNELS IN CELL ARRAY OR STRUCTURE,
%   PROVIDE CHANNEL NAME IN VARARGIN. EX: 'channel', 'ekar' *****
%
% OPTIONAL INPUTS
%   * lengthfilt - minimum trace length to consider for parameters
%                  DEFAULT = 0.
%
%   * twindow - restrict parameterization to a specified time window
%
%   * params - choose parameters from {'Mean';'Mean^2';'Max';'Mean Deriv';
%   'Frequency';'Mean Amp';'Mean Dur';'Max Amp';'IPI';'EndVal';'EndMax'}
%
%   * varnorm - normalize column by its own variance. DEFAULT = TRUE
%
% OUTPUTS
% parammat - cell x parameter matrix of parameter values
% 
%        p - contains parameter names (param), filtered tracks (badtrk),
%            parameter variance (varmap)
% ------------------------------------------------------------------------
% Version 2 

function [parammat, p] = ParameterizeData(d,z, varargin)
allParam = {'Mean';'Mean^2';'Max';'Mean Deriv';'Frequency';'Mean Amp';'Mean Dur';
    'Max Amp';'IPI';'EndVal';'Max^2';'EndMax'};
defaultParam = {'Mean';'Max';'Mean Deriv';'Frequency';'Mean Amp';'Mean Dur';
    'Max Amp';'IPI'};
% Set Default Parameters
p.params = defaultParam; p.lengthfilt = 0; p.varnorm = true; p.twindow = [];

% Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
% Splits pairs to a structure
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end

% parameter selection structure
fnames = strrep(allParam,'^','_'); fnames = strrep(fnames,' ','_'); % format names for struct
plog = ismember(lower(allParam),lower(p.params)); % logical for param selection
pulseparams = any(ismember(lower({'Frequency';'Mean Amp';'Mean Dur'; 'Max Amp';'IPI';}),...
    lower(p.params)));
for s = 1:numel(allParam) % pack into struct
   paramSelect.(fnames{s}) = plog(s); 
   IDX.(fnames{s}) = s;
end

%% * DATA FORMATING *
% If d is a cell array of traces turn it into a single channel
% NaN-padded matrix
if iscell(d) && isfield(d{1},'data') % ct_dataproc output format
    % check for old DEPRECATED dataproc format
    if size(d{1}.data,1) > 1
        error('Ew, old dataproc format, re-proc your data!')
    end
    % pack all data into one matrix
    d = d(~cellfun(@isempty,d)); % remove empty cells
    vcube = [];
    for i = 1:numel(d)
        vcube = [vcube; d{i}.data.(p.channel)]; % cat data from each XY
    end
    
elseif iscell(d) && size(d{1},2) == 1 % packed into structure w/o .data
    vcube = cell2mat(cellfun(@(x)x.(p.channel),d,'UniformOutput',false)');
    
else % Cell X Time matrix
    vcube = d;
end
% if z is a cell, make it a matrix of pulse data
if iscell(z)
    z = cat(1,z{1:end});
else
end

% Trim data matrices by time
if ~isempty(p.twindow)
    vcube = vcube(:,p.twindow);
else end

%% Calculate Parameters
% Initialize parameter matrix
parammat = nan(size(vcube,1),length(allParam));

% Initialize bad track matrix
p.badtrk = [];

% Loop through each row
for trk = 1:size(vcube,1)
    if sum(~isnan((vcube(trk,:)))) > p.lengthfilt
        trklength = sum(~isnan(vcube(trk,:)));
        
        if paramSelect.Mean % Mean of time window
            parammat(trk,IDX.Mean) = ...
                nanmean(nanmean(vcube(trk,:)));
        end
        
        if paramSelect.Mean_2 % Mean^2 (squared mean) 
            parammat(trk,IDX.Mean_2) = ...
                (nanmean(nanmean(vcube(trk,:))))^2;
        end
        
        if paramSelect.Max % Max of time window
            parammat(trk,IDX.Max) = ...
                max(max(vcube(trk,:)));
        end
        
        if paramSelect.Mean_Deriv % Mean Deriv
            tmp = diff(ct_filter(vcube(trk,:),'low',5));
            parammat(trk,IDX.Mean_Deriv) = ...
                nanmean(abs(tmp));
        end
        
        if paramSelect.EndVal% EndVal
            lst = find(~isnan(vcube(trk,:)),1,'last');
            parammat(trk,IDX.EndVal) = ...
                nanmean(vcube(trk,(lst-10):lst));
        end
        
        if paramSelect.Max_2  % Max^2
            parammat(trk,IDX.Max_2) = ...
                max(max(vcube(trk,:)))^2;
        end
        
        if paramSelect.EndMax % EndMax
            lst = find(~isnan(vcube(trk,:)),1,'last');
            parammat(trk,IDX.EndMax) = ...
                max(vcube(trk,(lst-10):lst));
        end
        % Pulse Parameters
        if pulseparams
            if ~isempty(z(trk).mpos) % Frequency
                if paramSelect.Frequency
                    parammat(trk,IDX.Frequency) = ...
                        numel(z(trk).mpos)/trklength;
                end
                if paramSelect.Mean_Amp % Mean Amplitude
                    parammat(trk,IDX.Mean_Amp) = ...
                        mean(z(trk).amp_mean);
                end
                if paramSelect.Mean_Dur % Mean Duration
                    parammat(trk,IDX.Mean_Dur) = ...
                        mean(z(trk).dur);
                end
                if paramSelect.Max_Amp% Max Peak Amplitude (Max Amp)
                    parammat(trk,IDX.Max_Amp) ...
                        = max(z(trk).amp_peak);
                end
                if paramSelect.IPI && numel(z(trk).mpos) > 2 % IPI
                    peakpos = z(trk).mpos;
                    peakdur = z(trk).dur;
                    parammat(trk,IDX.IPI) = ...
                        var((peakpos(2:end) + 0.5 * peakdur(2:end)) - ...
                        (peakpos(1:end-1) - 0.5 * peakdur(1:end-1)))/ ...
                        mean((peakpos(2:end) + 0.5 * peakdur(2:end)) - ...
                        (peakpos(1:end-1) - 0.5 * peakdur(1:end-1)));
                else parammat(trk,9) = 0;
                end
            elseif isempty(z(trk).mpos)
                varidx = ismember(allParam,{'Frequency','Mean Amp',...
                    'Mean Dur', 'Max Amp','IPI'});
                parammat(trk,varidx) = 0;
            end
        end
    else
        parammat(trk,1:numel(allParam)) = 0;
        p.badtrk = [p.badtrk; trk];
    end
end


parammat = parammat(:,struct2array(paramSelect));
% variance scale the data so all variables have a std of 1
if p.varnorm
    for i = 1:size(parammat,2)
        p.varmap(i) = std(parammat(:,i),'omitnan'); % divide by std
        if p.varmap(i) ~= 0
            parammat(:,i) = parammat(:,i)./p.varmap(i);
        else end
    end
else
end
end
