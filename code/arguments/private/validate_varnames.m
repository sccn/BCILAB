function validate_varnames(caller,varnames)
% Check for a given set of variable names for possible name clashes.
% validate_varnames(Caller,Varnames)
% 
% In:
%   Caller : name of the function (used in warning messages)
%   Varnames : cell array of variable names to check
% 
% This function checks for each variable name whether there is a function on the path with the same
% name; this is necessary because of a deficiency in how MATLAB handles variables that are assigned
% to a function's workspace using assignin (such variables are not first-class citizens in MATLAB).

    persistent already_checked;    
    if isempty(already_checked)
        already_checked = {}; end
    for n = fast_setdiff(varnames,already_checked)
        curname = n{1};
        try
            already_checked{end+1} = curname; %#ok<AGROW>
            existing_func = which(curname);
            if ~isempty(existing_func)
                if isempty(strfind(existing_func,'Java method'))
                    [path_part,file_part,ext_part] = fileparts(existing_func);
                    if ~any(strncmp('@',hlp_split(path_part,filesep),1))
                        % If this happens, it means that there is a function in one of the directories in
                        % MATLAB's path which has the same name as an argument of the specification. If this
                        % argument variable is copied into the function's workspace by arg_define, most MATLAB
                        % versions will (incorrectly) try to call that function instead of accessing the
                        % variable. I hope that they handle this issue at some point. One workaround is to use
                        % a longer argument name (that is less likely to clash) and, if it should still be
                        % usable for parameter passing, to retain the old name as a secondary or ternary
                        % argument name (using a cell array of names in arg()). The only really good
                        % solution at this point is to generally assign the output of arg_define to a
                        % struct.
                        disp([caller ': The argument name "' curname '" clashes with the function "' [file_part ext_part] '" in directory "' path_part '"; it is strongly recommended that you either rename the function or remove it from the path.']);
                    end
                else
                    % these Java methods are probably spurious "false positives" of the which() function
                    disp([caller ': There is a Java method named "' curname '" on your path; if you experience any name clash with it, please report this issue.']);
                end
            end
        catch e
            disp_once('Validation of the name %s failed; reason: %s',curname,e.message);
        end
    end
end
