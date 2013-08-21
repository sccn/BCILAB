function env_clear_classes
% Clear instances of BCI paradigm (and other) classes.
%
% This function is mainly useful when editing the code of BCI paradigms in BCILAB. Whenever after an
% edit a warning along the lines of "Cannot apply change to the class Paradigm*** because some
% instances of the old version still exist" comes up, this function may be called to clear all such
% stray instances without erasing the remaining state of BCILAB.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-12-29

global tracking;

% clear all micro-caches
hlp_microcache('clear');

% also clear all other stray class instances
try
    % make a backup of the tracking struct
    persistent backup; %#ok<TLEV>
    backup = tracking;
    % clear the classes, but keep this file (and the backup) locked
    mlock;
    now_clear;
    munlock;
    % restore the tracking variable from the backup
    restore_tracking(backup);
catch
    disp('Error while trying to clear classes from memory.');
end


function now_clear
% this needs to run in a different scope, otherwise the persistent variable ref would get lost
clear classes;


function restore_tracking(backup)
% this must be in a different scope, as we cannot have 2 global statements for the same variable
% in env_clear_memcaches()
global tracking;
tracking = backup;
