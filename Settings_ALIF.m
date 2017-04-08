function options = Settings_ALIF(varargin)
% 
% Settings_ALIF Constructs option structure for the algorithms
%
% EXAMPLES
%
% OPTIONS = SETTINGS_ALIF with no input arguments returns
%   setting structure with default values
%
% OPTIONS = SETTINGS_ALIF('NAME1',VALUE1,'NAME2',VALUE2,...) creates a
%   solution options structure OPTIONS in which the named properties have
%   the specified values.  Any unspecified properties have default values.
%   It is sufficient to type only the leading characters that uniquely
%   identify the property.  Case is ignored for property names.
%
% OPTIONS = SETTINGS_ALIF(OLDOPTIONS,'NAME1',VALUE1,...) alters an existing options
%   structure OLDOPTIONS.
%
%
% SETTINGS_ALIF PROPERTIES :
%
% GENERAL
% 
%   saveEnd          - (0) If value is 1 save outputs in VALUE.mat (only at termination)
%
%   verbose          - (1) 0 the method is silent, 1 normal, >1 loud, 
%
%   maxTime          - (Inf) Approximate max time in seconds one allows the
%                      algorithm to run. At the end of each iteration the  
%                      algorithm checks the elapsed time, if it has run for
%                      more than maxTime or maxTime - (time of last iteration)
%                      it stops iterating.
%                      Termination due to reaching of maxTime is signaled
%                      with a msg and in info.status. 
%       
%   plots            - (0) the algorithm does not produce plots,
%                      1 it produces them 
%
%   saveplots        - (0) not saving
%
% SPECIFIC PARAMETERS
% 
%  ALIF.delta     (0.0001)  
%  ALIF.ExtPoints (3)
%  ALIF.NIMFs     (100)
%  ALIF.xi        (1.6)
%
% ------------------------------------------------------
% EXAMPLE
%          
%   >> options = Settings_ALIF('ALIF.delta',0.08,'pl',1,'saveplots',1) 
%   >> IMF = ALIFv5(x,options)
%              
%  Executes algorithm ALIF with delta = 0.08, it produces plots and saves them                              
% ------------------------------------------------------      
%
% See also ALIFv5
%

% (Ripped from sdpsettings.m by Johan L�fberg)


% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    help Settings_ALIF
    return;
end


Names = {
    % General   
    'saveEnd'
    'verbose'
    'saveplots'
    'maxTime'
    'plots'
     % 'saveinIt'
     % 'logfile'
     
    
    % ALIF
    'ALIF.delta'
    'ALIF.ExtPoints'
    'ALIF.NIMFs'
    'ALIF.xi'
};

obsoletenames ={ % Use when options have become obsolete
};

[m,n] = size(Names);
names = lower(Names);

if (nargin>0) && isstruct(varargin{1})
    options = varargin{1};
    paramstart = 2;
else
    paramstart = 1;
    
    % General 
    options.saveinIt = 0;
    options.saveEnd = 0;
    options.verbose = 1;
    options.logfile = 1;
    options.maxTime = Inf;
    options.plots = 0;
    options.saveplots = 0; 
    
    % ALIF
    options.ALIF.delta = 0.0001; % used in the stopping criterion
    options.ALIF.ExtPoints=3;
    options.ALIF.NIMFs=100;
    options.ALIF.xi=1.6;
end

i = paramstart;
% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end
expectval = 0;       % start expecting a name, not a value
while i <= nargin
    arg = varargin{i};
    
    if ~expectval
        if ~ischar(arg)
            error(sprintf('Expected argument %d to be a string property name.', i));
        end
        
        lowArg = lower(arg);
        
        j_old = strmatch(lowArg,obsoletenames);
        if ~isempty(j_old)
            % For compability... No need yet
        end
        
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            error(sprintf('Unrecognized property name ''%s''.', arg));
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambiguous property name ''%s'' ', arg);
                msg = [msg '(' deblank(Names{j(1)})];
                for k = j(2:length(j))'
                    msg = [msg ', ' deblank(Names{k})];
                end
                msg = sprintf('%s).', msg);
                error(msg);
            end
        end
        expectval = 1;    % we expect a value next
    else
        eval(['options.' Names{j} '= arg;']);
        expectval = 0;
    end
    i = i + 1;
end

if expectval
    error(sprintf('Expected value for property ''%s''.', arg));
end

end

