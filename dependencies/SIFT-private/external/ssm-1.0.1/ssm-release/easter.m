function [m d] = easter(y)
% EASTER Easter date calculation for years 1583 to 4099
%
% y is a 4 digit year 1583 to 4099
%
% Using [m, d] = easter(y)
% d returns the day of the month of Easter
% m returns the month of Easter
%
% Using s = easter(y)
% s returns the Matlab serial date number for Easter
%
% Easter Sunday is the Sunday following the Paschal Full Moon
% (PFM) date for the year
%
% This algorithm is an arithmetic interpretation of the 3 step
% Easter Dating Method developed by Ron Mallen 1985, as a vast
% improvement on the method described in the Common Prayer Book
%
% Because this algorithm is a direct translation of the
% official tables, it can be easily proved to be 100% correct
%
% This algorithm derives values by sequential inter-dependent
% calculations, so ... DO NOT MODIFY THE ORDER OF CALCULATIONS!
%
% All variables are integer data types
%
% It's free!  Please do not modify code or comments!
% ==========================================================

%   Dim FirstDig, Remain19, temp    % intermediate results
%   Dim tA, tB, tC, tD, tE          % table A to E results

FirstDig    = fix(y/100);           % first 2 digits of year
Remain19    = mod(y, 19);           % remainder of year / 19

% calculate PFM date
temp        = fix((FirstDig - 15)/2) + 202 - 11*Remain19;
temp(FirstDig>26)   = temp(FirstDig>26) - 1;
temp(FirstDig>38)   = temp(FirstDig>38) - 1;
cond        = (FirstDig == 21) | (FirstDig == 24) | (FirstDig == 25) | ...
              (FirstDig == 33) | (FirstDig == 36) | (FirstDig == 37);
temp(cond)  = temp(cond) - 1;

temp        = mod(temp, 30);

tA          = temp + 21;
tA(temp==29)    = tA(temp==29) - 1;
tA(temp==28 & Remain19>10)  = tA(temp==28 & Remain19>10) - 1;

% find the next Sunday
tB          = mod(tA - 19, 7);

tC          = mod(40 - FirstDig, 4);
tC(tC==3)   = tC(tC==3) + 1;
tC(tC>1)    = tC(tC>1) + 1;

temp        = mod(y, 100);
tD          = mod(temp + fix(temp/4), 7);

tE          = mod(20 - tB - tC - tD, 7) + 1;
d           = tA + tE;

% return the date
m           = repmat(3, size(d));
m(d>31)     = 4;
d(d>31)     = d(d>31) - 31;
