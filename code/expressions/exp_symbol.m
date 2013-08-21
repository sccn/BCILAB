function func = exp_symbol(name)
% Declare a symbol with the given name.
% Symbol = exp_symbol(Name)
% 
% In: 
%   Name   : name of the symbol to be declared; must be a valid MATLAB identifier (a letter followed
%            by alphanumeric chars and/or underscores)
%
% Out:
%   Symbol : a symbol that can be used inside expressions
%
% Examples:
%   f = exp_symbol('f')
%   f(1,2,f(4,5,6),3)  --> gives an expression
%   
%   sin = exp_symbol('sin')
%   sin(0.5) --> gives an expression that can be evaluated into a number using exp_eval()
%
% See also:
%   exp_eval, exp_declare_symbols
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-19

persistent symboltable; % a cache of recently produced symbols... (these are expensive to create!)

if ~exist('name','var') || ~ischar(name)
    error('A name must be specified to identify the symbol.'); end

if isfield(symboltable,name)
    % take the declaration from the symbol table
    func = symboltable.(name);
else    
    % evaluating the function declaration with name substituted literally (rather than as a variable
    % reference) keeps the name visible/accessible from the literal char() transcription of the
    % lambda function
    func = eval(['@(varargin)struct(''head'',{@' name '},''parts'',{varargin})']);
    symboltable.(name) = func;
end
