% HandleVar class
% A simple, generic handle class that can be used to create shared
% variables. All copies of a HandleVar variable reference the same
% underlying object. Only has a single property ('data'), which can be
% accessed freely.

classdef HandleVar < handle
    properties
        data;
    end
end