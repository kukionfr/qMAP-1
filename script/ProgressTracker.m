classdef ProgressTracker < handle
    properties
        Value = 0;
        Total = 1;
        WaitbarHandle
    end
    methods
        function obj = ProgressTracker(total, h)
            obj.Total = total;
            obj.WaitbarHandle = h;
        end
        function increment(obj)
            obj.Value = obj.Value + 1;
            waitbar(obj.Value / obj.Total, obj.WaitbarHandle, ...
                sprintf('Processing pair %d of %d...', obj.Value, obj.Total));
        end
    end
end
