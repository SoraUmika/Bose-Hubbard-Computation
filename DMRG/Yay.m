classdef Yay < handle
    properties
        m
        refcnt;
    end
    
    methods
        function self = Yay(m)
            self.m = m;
            self.refcnt = 1;
        end
    end
end

