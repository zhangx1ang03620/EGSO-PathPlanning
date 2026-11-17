classdef FobjCounter < handle
    % FobjCounter  Wraps an objective function and counts calls, for NFE
    % (number of function evaluations) reporting. Note: calls made
    % directly to fitness() (e.g. inside LocalSearch.m) are not seen.
    properties
        count = 0
        inner
    end
    methods
        function obj = FobjCounter(fobj)
            obj.inner = fobj;
        end
        function f = call(obj, x)
            obj.count = obj.count + 1;
            f = obj.inner(x);
        end
    end
end
