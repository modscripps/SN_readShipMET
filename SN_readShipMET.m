function MET = SN_readShipMET(filename)
% SN_READSHIPMET - reads SIO ships MET data files
%
% SN_READSHIPMET(FILENAME) returns a MET structure of variables described
% for MET data files
%
% Created 2012/03/30 San Nguyen
% Updated 2012/06/09 San Nguyen - updated to read MET data for most ships
% not just the ones on the R/V Revelle
% Updated 2015/09/02 San Nguyen - update to combine MET data without losing
% any fields added later on
%

if ischar(filename)
    switch exist(filename)
        case 2 % if it is a file
            fid = fopen(filename,'r');
            MET = SN_readShipMET(fid);
            fclose(fid);
        case 7
            my_METfiles = dir([filename '/*.MET']);
            if isempty(my_METfiles)
                MET = [];
                return
            else
                disp(['reading ' my_METfiles(1).name]);
                MET = SN_readShipMET([filename '/' my_METfiles(1).name]);
                varnames = fieldnames(MET);
                MET = repmat(MET,size(my_METfiles));
                for i = 2:length(my_METfiles)
                    disp(['reading ' my_METfiles(i).name]);
                    MET_tmp = SN_readShipMET([filename '/' my_METfiles(i).name]);
                    for k = 1:numel(varnames)
                        if strcmpi(varnames{k},'README')
                            MET(i).(varnames{k}) = MET_tmp.(varnames{k});
                        end
                        if isfield(MET_tmp,varnames{k})
                            MET(i).(varnames{k}) = MET_tmp.(varnames{k});
                        else
                            MET(i).(varnames{k}) = NaN(size(MET_tmp.(varnames{1})));
                        end
                    end
                end
                MET = SN_combineMET(MET);
            end
        otherwise
            error('MATLAB:SN_readShipMET:wrongFILENAME','Invalid file specified.');
    end
    
elseif iscellstr(filename)
    MET = SN_readShipMET(filename{1});
    MET = repmat(MET,length(filename),1);
    for i = 2:length(filename)
        MET(i) = SN_readShipMET(filename{i});
    end
    MET = SN_combineMET(MET);
else
    
    if (filename<1)
        error('MATLAB:SN_readShipMET:wrongFID','FID is invalid');
    end
    
    MET = SN_readShipMETFile(filename);
    
    return
end

end

% reading MET files through FID
function MET = SN_readShipMETFile(fid)
% struct_vars = {...
%     'Time'; 'AT'; 'BP'; 'BC'; 'BS'; 'RH'; 'RT';...
%     'DP'; 'PR'; 'LB'; 'LT'; 'LW'; 'SW'; 'PA'; 'WS';...
%     'WD'; 'TW'; 'TI'; 'WS_2'; 'WD_2'; 'TW_2'; 'TI_2';...
%     'TT'; 'TC'; 'SA'; 'SD'; 'SV'; 'TG'; 'FI'; 'PS';...
%     'TT_2'; 'TC_2'; 'SA_2'; 'SD_2'; 'SV_2'; 'TG_2';...
%     'OC'; 'OT'; 'OX'; 'OS'; 'FL'; 'FI_2'; 'PS_2';...
%     'VP'; 'VR'; 'VH'; 'VX'; 'VY'; 'GY'; 'MB'; 'BT';...
%     'SH'; 'SM'; 'SR'; 'LA'; 'LO'; 'GT'; 'CR'; 'SP';...
%     'ZD'; 'GA'; 'GS'; 'ZO'; 'ZS'; 'ZT'; 'XX'; 'ZO_2';...
%     'ZS_2'; 'ZT_2'; 'XX_2'; 'PZ'; 'PZ_2'; 'IP'; 'IV';...
%     'IA'; 'SW_2'; 'SW_3'};

frewind(fid)
fgetl(fid);
str = fgetl(fid);
METdate = floor(datenum(str(7:end),'dd-mmm-yy  HH:MM:SS'));

fgetl(fid);
str = fgetl(fid);

str = strrep(str,'#','');
str = strrep(str,'-','_');
struct_vars = textscan(str,'%s','EndOfLine',' ');
struct_vars = struct_vars{1};


data_str = char(ones(1,2*length(struct_vars)));
time_ind = NaN;
for i = 1:length(struct_vars)
    if strcmpi(struct_vars{i},'Time')
        data_str((i-1)*2+(1:2)) = '%s';
        time_ind = i;
    else
        data_str((i-1)*2+(1:2)) = '%f';
    end
end

if isnan(time_ind)
    error('MATLAB:SN_readShipMET','We can not find the time in the data file. Please check');
end

frewind(fid)
data = textscan(fid,[data_str '%*[^\n]'],'delimiter',' \t','MultipleDelimsAsOne',true,'CommentStyle','#');
data{time_ind} = datenum(data{time_ind},'HHMMSS')-datenum('000000','HHMMSS')+METdate;

%make sure that the time from one day is not mistaken for the other day
time_inv_ind = find(diff(data{time_ind})<0);
data{time_ind}(time_inv_ind) = data{time_ind}(time_inv_ind) - 1;

MET = struct();

for i = 1:length(struct_vars)
    MET.(struct_vars{i}) = data{i};
end

MET.README = {...
    'AT:    Air Temperature';
    'RH:    Relative Humidity';
    'RT:    Relative Temperature';
    'DP:    Dew Point';
    'BP:    Barometric Pressure';
    'BC:    Barometric Compression Temperature';
    'PR:    Precipitation';
    'PT:    Precipitation Temperature';
    'WS:    Wind Speed';
    'WD:    Wind Direction';
    'WK:    ';
    'TK:    ';
    'TW:    ';
    'TI:    ';
    'SW:    Short wave radiation';
    'LW:    Long wave radiation';
    'LD:    Long wave radation dome temperature';
    'LB:    Long wave radation body temperature';
    'LT:    Long wave radation thermopile voltage';
    'FM:    Flowmeter';
    'FI:    Flowmeter';
    'ST:    Surface seawater temperature';
    'SC:    Surface seawater conductivity';
    'VT:    A/D Volts';
    'FL:    Fluorometer';
    'TB:    ';
    'TR:    Transmissometer';
    'BA:    ';
    'OC:    Oxygen';
    'OT:    ';
    'OS:    ';
    'OX:    ';
    'OG:    ';
    'TT:    Thermosalinograph';
    'TC:    ';
    'SA:    ';
    'SD:    ';
    'SV:    ';
    'XX:    ';
    'PA:    Surface PAR';
    'WT:    Auxiliary water temperature';
    'AX:    Auxiliary air temperature';
    'IP:    Instrumentation';
    'IT:    ';
    'IS:    ';
    'IA:    ';
    'IV:    ';
    'IX:    ';
    'PS:    Pressure';
    'BT:    Water depth';
    'LA:    NMEA messages';
    'LO:    ';
    'GT:    ';
    'CR:    ';
    'SP:    ';
    'ZD:    ';
    'GA:    ';
    'GS:    ';
    'GY:    Gyro';
    'SH:    Ashtech heading';
    'SM:    Ashtech pitch';
    'SR:    Ashtech roll';
    'TS:    Time server';
    'ZO:    Winch wire out';
    'ZS:    Winch wire speed';
    'ZT:    Wihch wire tension';
    'PH:    Alkalinity (pH)';
    'VP:    Vertical reference unit pitch';
    'VR:    Vertical reference unit roll';
    'VH:    Vertical reference unit heave';
    'VY:    Vertical reference unit ';
    'VX:    Vertical reference unit ';
    'SL:    Ship''s speed log (speed over water)';
    'PZ:    Uncontaiminated seawater pump status';
    };
end

function MET = SN_combineMET(varargin)
% SN_combineMET - combines MET data files in MATLAB format that was
% converted using SN_readShipMET
%
% SN_combineMET(MET1,MET2,MET3,...) returns a MET structure of variables described
% for MET data files
%
% Written 2015/09/02 - San Nguyen snguyen@opg1.ucsd.edu

if nargin < 1 
    MET = [];
    return;
end

if nargin == 1
   MET = varargin{1};
   if length(MET) == 1
       return;
   end
   evalstr = 'MET = SN_combineMET(';
   for i=1:(length(MET)-1)
       evalstr = strcat(evalstr, 'MET(', num2str(i), '),');
   end
   evalstr = strcat(evalstr, 'MET(', num2str(length(MET)), '));');
   eval(evalstr);
   return
end

MET_fields = fieldnames(varargin{1});
for i = 2:nargin
    tmp_fields = fieldnames(varargin{i});
    for j = 1:numel(tmp_fields)
        if ~ismember(tmp_fields{j},MET_fields)
            MET_fields{end+1} = tmp_fields{j};
        end
    end
end

for i=1:(length(MET_fields))
    if strcmpi(MET_fields{i},'README')
        continue;
    end
    evalstr = strcat('MET.', MET_fields{i}, '= [');
    for j=1:(nargin-1)
        if ~isfield(varargin{j},(MET_fields{i}))
            varargin{j}.(MET_fields{i}) = NaN(size(varargin{j}.Time));
        end
        evalstr = strcat(evalstr, 'varargin{', num2str(j), '}.', MET_fields{i}, ';');
    end
    if ~isfield(varargin{nargin},(MET_fields{i}))
        varargin{nargin}.(MET_fields{i}) = NaN(size(varargin{nargin}.Time));
    end
    evalstr = strcat(evalstr, 'varargin{', num2str(nargin), '}.', MET_fields{i}, '];');
    eval(evalstr);
end

if (isfield(varargin{1},'README'))
    MET.README = varargin{1}.README;
end

% sort out time
[~, I] = sort(MET.Time);

for i=1:(length(MET_fields))
    if strcmpi(MET_fields{i},'README')
        continue;
    end
    MET.(MET_fields{i}) = MET.(MET_fields{i})(I);
end

%find unique time
[~, I, ~] = unique(MET.Time);

for i=1:(length(MET_fields))
    if strcmpi(MET_fields{i},'README')
        continue;
    end
    MET.(MET_fields{i}) = MET.(MET_fields{i})(I);
end

end