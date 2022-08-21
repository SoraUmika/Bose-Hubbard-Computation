classdef InnerObject < handle   
    properties
        site_indices;
        operator_names;
        area
        ID
    end
    methods
        function self = InnerObject(site_indices, operator_names, area)
            self.site_indices = site_indices;
            self.operator_names = operator_names;
            self.area = area;
            self.ID = strcat(mat2str(site_indices),string(operator_names));
        end
    end
end
