% parse_args()
% Parse user arguments supplied as name/value pairs. Return arguments as a
% struct, setting default values where not specified by user.
%
% Usage:
%   [ args ] = parse_args( user_args, default_args, match_case )
%
% Inputs:
%   user_args
%     A cell array containing name/value pairs. Typically this is passed to
%     a function as varargin.
%   default_args
%     A struct containing a field for each argument. Corresponding values
%     are the default values.
%   match_case (Optional)
%     If true, argument names are case sensitive. [Default: false]
%
% Outputs:
%   args
%     A struct containing a field for every argument. Corresponding values
%     are taken from user_args if present, otherwise default_args.
%     Case of field names matches default_args.
%
% Notes:
%   Throws an error if user-supplied arguments aren't valid name/value
%   pairs or contain an argument not in default_args.

function [ args ] = parse_args( user_args, default_args, match_case )

    % case insensitive by default
    if nargin < 3
        match_case = false;
    end

    % parse user args into name/value pairs
    user_names = user_args(1:2:end);
    user_vals = user_args(2:2:end);
    if ~iscellstr(user_names) || numel(user_names) ~= numel(user_vals)
        error('Arguments must be valid name/value pairs');
    end
    
    default_names = fieldnames(default_args);
    
    % set args to defaults, replace w/ user-specified values where given
    args = default_args;
    unknown_names = {};
    for it = 1 : numel(user_names)
        if match_case
            idx = find(strcmp(user_names{it}, default_names));
        else
            idx = find(strcmpi(user_names{it}, default_names));
        end
        if isempty(idx)
            unknown_names{end+1} = user_names{it};
        else
            args.(default_names{idx}) = user_vals{it};
        end
    end
    
    % make sure user didn't specify unknown args
    if ~isempty(unknown_names)
        error(sprintf( ...
            'Unknown arguments: %s', ...
            strjoin(unknown_names, ', ') ...
        ));
    end
%
end