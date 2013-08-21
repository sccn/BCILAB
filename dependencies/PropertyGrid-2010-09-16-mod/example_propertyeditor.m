% Demonstrates how to use the property editor.
%
% See also: PropertyEditor

% Copyright 2010 Levente Hunyadi
function example_propertyeditor

% create figure
f = figure( ...
    'MenuBar', 'none', ...
    'Name', 'Property editor demo - Copyright 2010 Levente Hunyadi', ...
    'NumberTitle', 'off', ...
    'Toolbar', 'none');
items = { SampleObject SampleObject };
editor = PropertyEditor(f, 'Items', items);
editor.AddItem(SampleNestedObject, 1);
editor.RemoveItem(1);
