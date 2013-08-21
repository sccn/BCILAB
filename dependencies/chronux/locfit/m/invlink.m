function y0 = invlink(y,fali)
% inverse link function for locfit.
% y is a vector of raw fitted values.
% fali is the integer [family link] vector from locfit.
% output is the inv. link.

link = fali(2);

switch(link)
    case 3   % identity
        y0 = y;
    case 4   % log
        y0 = exp(y);
    case 5   % logit - should invert carefully!
        y0 = 1 - 1./(1+exp(y));
    case 6   % inverse
        y0 = 1/y;
    case 7   % sqrt
        y0 = y*abs(y);
    case 8   % arcsin
        y0 = sin(y)*sin(y);
    otherwise
        disp('invlink: Unknown link function');
end;
