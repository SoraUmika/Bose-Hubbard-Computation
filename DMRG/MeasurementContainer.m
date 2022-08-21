classdef MeasurementContainer < handle
    properties
        measurements
        sites_by_area
        site_class
        measurements_by_site
        processed_measurements = MeasurementWrapper.empty
        g_some_dict
        trmat_disk
        DMRG_obj
        returned_measurements;
    end
    
    methods
        function self = MeasurementContainer(DMRG_obj, measurements)
            self.DMRG_obj = DMRG_obj;          
            self.measurements = measurements;
            
            %We first order the measurements by section, or area
            L = self.DMRG_obj.model.L;
            [LEFT_BLOCK, LEFT_SITE, RIGHT_BLOCK, RIGHT_SITE] = deal(1,2,3,4);
            self.sites_by_area = containers.Map('KeyType','uint64','ValueType','any');
            if self.DMRG_obj.model.bound_cond == BoundCond.open
                self.sites_by_area(LEFT_BLOCK) = 1:(floor(L/2)-1);
                self.sites_by_area(LEFT_SITE) = floor(L/2);
                self.sites_by_area(RIGHT_BLOCK) = flip((floor(L/2)+2):L);
                self.sites_by_area(RIGHT_SITE) = floor(L/2)+1;
            elseif DMRG_obj.model.bound_cond == BoundCond.periodic
                self.sites_by_area(LEFT_BLOCK) = 1:(floor(L/2)-1);
                self.sites_by_area(LEFT_SITE) = floor(L/2);
                self.sites_by_area(RIGHT_BLOCK) = flip(floor(L/2+1):(L-1));
                self.sites_by_area(RIGHT_SITE) = L;                
            end
            self.site_class = containers.Map('KeyType','uint64','ValueType','uint64');
            areas = cell2mat(keys(self.sites_by_area));
            for i=1:self.sites_by_area.Count
                area_indices = self.sites_by_area(areas(i));
                for j=1: length(area_indices)
                    self.site_class(area_indices(j)) = i;
                end
            end
            
            %Now we order the measurements by site
            self.measurements_by_site = containers.Map('KeyType','uint64','ValueType','any');
            for i=1:length(measurements)
                meas_desc = measurements(i);
                operator_count = length(meas_desc.site_indices);
                [meas_desc.site_indices, sorted_indices] = sort(meas_desc.site_indices);
                meas_desc.operator_names = meas_desc.operator_names(sorted_indices);

                site_indices = [];
                measWrapper = MeasurementWrapper();
                for j=1:operator_count
                    site_index = meas_desc.site_indices(j);                 
                    measWrapper.measurement(j) = OuterObject(site_index,meas_desc.operator_names(j),self.site_class(site_index),false);
                    site_indices(j) = site_index;
                end
                self.processed_measurements(i) = measWrapper;
                site_indices = unique(site_indices);
                for j=1:length(site_indices)
                    site_index = site_indices(j);
                    measurements_on_site = measWrapper;
                    if isKey(self.measurements_by_site, site_index)
                        measurements_on_site = self.measurements_by_site(site_index);
                        measurements_on_site(end+1) = measWrapper;
                    end
                    self.measurements_by_site(site_index) = measurements_on_site;
                end              
            end
            assert(length(measurements) == length(self.processed_measurements));
            
            map1 = containers.Map('KeyType','char','ValueType','any');
            map2 = containers.Map('KeyType','char','ValueType','any');
            map3 = containers.Map('KeyType','char','ValueType','any');
            map4 = containers.Map('KeyType','char','ValueType','any');
            self.g_some_dict = {map1, map2, map3, map4};
            
            lb_msize = self.BuildBlock(LEFT_BLOCK, 'l');
            rb_msize = self.BuildBlock(RIGHT_BLOCK, "r");
            
            LEFT_SITE_AREA = self.sites_by_area(LEFT_SITE);
            self.HandleOperatorsOnSite(LEFT_SITE_AREA, LEFT_SITE, self.g_some_dict{LEFT_SITE}, 1);

            RIGHT_SITE_AREA = self.sites_by_area(RIGHT_SITE);
            self.HandleOperatorsOnSite(RIGHT_SITE_AREA, RIGHT_SITE, self.g_some_dict{RIGHT_SITE}, 1);
            
            rpsi0 = self.DMRG_obj.gnd_state(self.DMRG_obj.rbi);
            
            d = self.DMRG_obj.model.d;
            mm_orig = {speye(lb_msize), speye(d), speye(rb_msize),speye(d)};
            self.returned_measurements = zeros(length(self.processed_measurements),1);
            
            fprintf("Taking measurements...\n");
            for i=1: length(self.processed_measurements)
                pm = self.processed_measurements(i).measurement;
                mm = mm_orig;
                for j=1: length(pm)
                    obj = pm(j);
                    some_dict = self.g_some_dict{obj.obj.area};
                    yay_obj = some_dict(obj.obj.ID);
                    mm{obj.obj.area} = yay_obj.m;
                end
                big_m = kron(kron(mm{1},mm{2}), kron(mm{3},mm{4}));
                big_m = big_m(self.DMRG_obj.rbi, self.DMRG_obj.rbi);
                
                ev = (ctranspose(rpsi0))*(big_m*rpsi0);
                fprintf("Measurement %d: %d\n", i, ev(1));
                self.returned_measurements(i) = ev(1);
            end
        end
        
        function HandleOperatorsOnSite(self, site_index, area, some_dict, msize)
            measuremets_in_site = self.measurements_by_site(site_index);
            for measurement_num=1:length(measuremets_in_site)
                measWrapper = measuremets_in_site(measurement_num);
                measurement = measWrapper.measurement;
                for i=1: length(measurement)
                    obj = measurement(i);
                    if obj.obj.site_indices(1) ~= site_index
                        continue;
                    end
                    assert(obj.built == false)
                    if isKey(some_dict,obj.obj.ID)
                        yay_obj = some_dict(obj.obj.ID);
                        yay_obj.refcnt = yay_obj.refcnt + 1;
                        %some_dict(ID) = yay_obj;
                    else
                        assert(length(obj.obj.operator_names) == 1)
                        sso = self.DMRG_obj.model.sso(obj.obj.operator_names(1));
                        mat = kron(speye(msize),sso);
                        some_dict(obj.obj.ID) = Yay(mat);
                    end          
                    measurement(i).built = true;
                end 
                        
                for i=(length(measurement)-1):-1:1
                    areaCond1 = measurement(i).obj.area == area;
                    areaCond2 = measurement(i+1).obj.area == area;
                    areaCond = areaCond1 && areaCond2; 
                    if measurement(i).built && measurement(i+1).built && areaCond
                        c1 = measurement(i);
                        c2 = measurement(i+1);
                                              
                        o1 = some_dict(c1.obj.ID);
                        o1.refcnt = o1.refcnt - 1;
                        some_dict(c1.obj.ID) = o1;
                        if o1.refcnt == 0
                            remove(some_dict, c1.obj.ID);
                        end
                  
                        o2 = some_dict(c2.obj.ID);
                        o2.refcnt = o2.refcnt - 1;
                        some_dict(c2.obj.ID) = o2;
                        if o2.refcnt == 0
                            remove(some_dict, c2.obj.ID);
                        end
                        c3 = OuterObject([c1.obj.site_indices, c2.obj.site_indices], ...
                                         [c1.obj.operator_names, c2.obj.operator_names], ...
                                         area , true);
                           
                        if isKey(some_dict, c3.obj.ID)
                            yay_obj = some_dict(c3.obj.ID);
                            yay_obj.refcnt = yay_obj.refcnt + 1;                        
                        else
                            yay_obj = Yay(o1.m*o2.m);
                        end
                        some_dict(c3.obj.ID) = yay_obj;
                        
                        measurement(i+1) = [];
                        measurement(i) = c3;
                        measWrapper.measurement = measurement;
                    end
                end
            end
        end
        
        function msize = BuildBlock(self, area, block_label)
            assert(block_label == 'r' || block_label == 'l');
            some_dict = self.g_some_dict{area};
            assert(some_dict.Count == 0);
            msize = 1;
            
            sites_in_area = self.sites_by_area(area);
            for i=1: length(sites_in_area)
                site_index = sites_in_area(i);
                some_dict_keys = keys(some_dict);
                for j=1: length(some_dict_keys)
                    k = some_dict_keys(j);
                    k = k{1};
                    v = some_dict(k);
                    v.m = kron(v.m, speye(self.DMRG_obj.model.d));
                end
                
                self.HandleOperatorsOnSite(site_index, area, some_dict, msize);  
                
                if i==1
                    msize = self.DMRG_obj.model.d;
                else
                    if block_label == 'r'
                        T_matrix = self.DMRG_obj.r_T_matrices(i);
                    else
                        T_matrix = self.DMRG_obj.l_T_matrices(i);
                    end
                    some_dict_keys = keys(some_dict);
                    for j=1: length(some_dict_keys)
                        k = some_dict_keys(j);
                        k = k{1};
                        v = some_dict(k);
                        assert(v.refcnt > 0);
                        %assert k.area == area
                       
                        v.m = sparse(T_matrix'*v.m*T_matrix);
                    end
                    [row, msize] = size(T_matrix);
                end
            end
            
        end
        
    end
end
