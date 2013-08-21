function output = resetintegrator(input, reset)
%RESETINTEGRATOR   Running counter with resets.
%   OUTPUT = RESETINTEGRATOR(INPUT, RESET), where INPUT and RESET are 
%   length N vectors, returns a length N vector OUTPUT equal to the
%   cumulative sum of INPUT, reset to zero on each sample k for which
%   RESET(k) == 0.   
%
%   OUTPUT = RESETINTEGRATOR(INPUT) uses INPUT as the RESET vector; it
%   is equivalent to OUTPUT = RESETINTEGRATOR(INPUT, INPUT).

if (nargin < 2),  reset = input;  
elseif (length(input) ~= length(reset)), 
	error('INPUT and RESET vectors must be of the same length.');
end;
if (isempty(input)), output = zeros(size(input));  return;  end;  % trap special case

output = CORE_resetintegrator(double(input), logical(reset));
