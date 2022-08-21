classdef OuterObject < handle
    properties
        obj
        built
    end
    
    methods
        function self = OuterObject(site_indices, operator_names, area, built)
           self.obj = InnerObject(site_indices, operator_names, area);
           self.built = built;
        end
    end
end
