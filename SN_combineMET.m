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