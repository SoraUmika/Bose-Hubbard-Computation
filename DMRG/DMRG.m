classdef DMRG < handle
    properties
        model
        l_block_list
        l_T_matrices
        r_T_matrices
        r_block_list
        
        gnd_state 
        energy
        rbi     
        
        logging = false;
    end
    
    methods
        function self = DMRG(model)            
            self.model = model;            
            self.l_block_list = containers.Map('KeyType','uint64','ValueType','any');
            self.l_T_matrices = containers.Map('KeyType','uint64','ValueType','any');
            self.r_block_list = containers.Map('KeyType','uint64','ValueType','any');
            self.r_T_matrices = containers.Map('KeyType','uint64','ValueType','any');          
        end
%--------------------------------------------------------------------------        
        %this function performs a single DMRG steps that involves:
        %enlarging the sys and env block; constructing the superblock
        %hamiltonian; get the reduce basis indices using conserved quantum
        %numbers or symmetries; finding the reduce-density-matrix of the 
        %sys blockand its corresponding eigen values and vectors; using 
        %those eigenvalue to construct a set of new basis for truncation;
        %updates the the operators associated with the sys block.
        function [sys_block, gnd_state, energy, rbi, T_matrix] = ...
                SingleDMRGStep(self, sys_block, env_block, m, direction, conserved_QNum)    
            sys_block = self.model.EnlargeBlock(sys_block, direction);          
            env_block = self.model.EnlargeBlock(env_block, direction);                      
            superBlock_H = self.model.ConstructSuperBlock_H(sys_block, env_block);
            
        %------------------------------------------------------------------     
            %basis_QNum is an arary that contains the value of all possible
            %quantum number for that block. The following KeyValueMap maps
            %the quantum numbers to the indices that they were found. For
            %example, consider the sys_QNum_Indices_Map. If the quantum number 
            %3 is found at indices 3,7,2 then sys_QNum_Indices_Map(3)
            %will give us {[3,7,2]}, note that this is in a cell array. 
            sys_QNum_Indices_Map = KeyValueMap(sys_block.basis_QNum);
            env_QNum_Indices_Map = KeyValueMap(env_block.basis_QNum);
            
            %Now we consider the restriction of conserved quantum numbers
            sys_possible_QNum_List = cell2mat(keys(sys_QNum_Indices_Map));
            sys_QNum_Count = sys_QNum_Indices_Map.Count;
            new_sys_QNum_Indices_Map = ... 
                containers.Map(sys_possible_QNum_List, cell(1,sys_QNum_Count));             
            
            rbi = [];  
            rbi_index = 1;
            for i=1: sys_QNum_Count
                sys_QNum = sys_possible_QNum_List(i);
                %here we calculate env_QNum by assuming the QNum on the sys_block
                env_QNum = conserved_QNum - sys_QNum;
                %checks whether the calculated env_QNum is in the possible
                %QNum of the env_block. 
                if isKey(env_QNum_Indices_Map, env_QNum)
                    %gets the indices(basis) that matches whatever sys_QNum is
                    sys_QNum_Indices = sys_QNum_Indices_Map(sys_QNum);
                    for j=1: length(sys_QNum_Indices)
                        index_offset = env_block.basis_size*(sys_QNum_Indices(j)-1);
                        %gets the indices(basis) that matches whatever env_QNum is
                        env_QNum_Indices = env_QNum_Indices_Map(env_QNum);
                        for k=1: length(env_QNum_Indices)
                            rbi(rbi_index) = index_offset + env_QNum_Indices(k);                                       
                            tmp_arr = new_sys_QNum_Indices_Map(sys_QNum);
                            new_sys_QNum_Indices_Map(sys_QNum) = [tmp_arr rbi_index];
                            
                            rbi_index = rbi_index + 1;
                        end
                    end
                end
            end                         
            restricted_superBlock_H = ... 
                superBlock_H(rbi,rbi);
            [restricted_gnd_state, energy] = eigs(restricted_superBlock_H,1,'smallestreal');
            
            gnd_state = zeros(sys_block.basis_size*env_block.basis_size,1);     
            for i=1: length(rbi)
               z = rbi(i);
               gnd_state(z) = restricted_gnd_state(i);
            end
            
        %------------------------------------------------------------------       
            %we now calculate the reduce density matrix. Here, we make them in
            %block matrix form where each block of the matrix represents the 
            %corresponding quantum numbers 
            %rho_e_values = zeros(1,sys_block.basis_size);        
            new_sys_QNum = keys(new_sys_QNum_Indices_Map);
            new_sys_Indices = values(new_sys_QNum_Indices_Map);   

            rho_e_values = [];
            rho_e_vects = zeros(sys_block.basis_size,1);   
            new_basis_QNum = [];
            current_evalue_index = 1;
            for i=1:new_sys_QNum_Indices_Map.Count
                if ~isempty(new_sys_Indices{i})
                    psi0_QNum = restricted_gnd_state(new_sys_Indices{i});          
                    sys_QNum_Indices = sys_QNum_Indices_Map(new_sys_QNum{i}).'; 
                    psi0_QNum = reshape(psi0_QNum, [],length(sys_QNum_Indices)).';      
                    rho_block = psi0_QNum*psi0_QNum';           
                    [e_vects, diag_mat] = eig(rho_block);  
                    e_values = diag(diag_mat);
                    current_QNum_Indices = sys_QNum_Indices_Map(new_sys_QNum{i});
                    for j = 1: length(e_values)
                        rho_e_values(current_evalue_index) = e_values(j);
                        rho_e_vects(current_QNum_Indices,current_evalue_index) = e_vects(:,j);
                        new_basis_QNum(current_evalue_index) = new_sys_QNum{i};
                        current_evalue_index = current_evalue_index + 1;
                    end
                end
            end
        %------------------------------------------------------------------  
            %reorder the eigralues, eigen vectors, the corresponding quantum 
            %number in descending eigenvalue order. 
            [rho_e_values, Index] = sort(rho_e_values, 'descend');
            rho_e_vects = rho_e_vects(:,Index);        
            new_basis_QNum = new_basis_QNum(Index);

            %truncation 
            NKeep = min(length(rho_e_values), m);
            T_matrix = zeros(sys_block.basis_size, NKeep);
            T_matrix = rho_e_vects(:, 1:NKeep);
            new_basis_QNum = new_basis_QNum(:, 1:NKeep);
            trunc_err = 1 - sum(rho_e_values(1:NKeep)); 
        %------------------------------------------------------------------              
            %print out the chain length, energy and the associated
            %truncation error associated with this iteration
            total_L = sys_block.length+env_block.length;
            if self.logging
                fprintf("L=%i, E=%d, err=%d\n", total_L, energy, trunc_err); 
            end           
        %------------------------------------------------------------------              
            %update our curent block
            sys_block.basis_size = NKeep;
            sys_block.length = sys_block.length;           
            sys_block.basis_QNum = new_basis_QNum;
            
            new_op_name = keys(sys_block.op_list); 
            for i=1: length(new_op_name)
               op_name = new_op_name(i);
               op_mat = sys_block.op_list(op_name{1});
               sys_block.op_list(op_name{1}) = sparse(T_matrix'*op_mat*T_matrix);
            end             
        end

%--------------------------------------------------------------------------
        %performs iDMRG
        function [energy, gnd_state] = iDMRG(self, dmrg_input)
            assert(isfield(dmrg_input, "m_warmup"), "please specify m_warmup in dmrg_input");
            fprintf("------------------------------------------------- \n");
            chain_length = self.model.L;
            m = dmrg_input.m_warmup;
            conserved_QNum = self.model.conserved_QNum;
            
            l_block = self.model.InitSingleSiteBlock();   
            iteration_count = floor((chain_length-2)/2);
            direction = 'r';
            for iteration=1: iteration_count 
                r_block = l_block;
                [l_block, gnd_state, energy, ~, ~] = ...
               self.SingleDMRGStep(l_block, r_block, m, direction, conserved_QNum);
            end
            
            %checks whether our chain length is odd or even, if its odd,
            %then our new left block will be 1 sites larger than the right
            %block
            if mod(chain_length,2) ~= 0
                [~, gnd_state, energy, ~, ~] = ...
               self.SingleDMRGStep(l_block, r_block, m, direction, conserved_QNum);             
            end
            
            fprintf("infinite DRMG success with L=%d, E=%d \n", chain_length, energy);
            fprintf("------------------------------------------------- \n");
        end
%--------------------------------------------------------------------------        
        %performs fDMRG by first building up the chain using iDMRG, then
        %the sweep.
        function [energy, gnd_state] = fDMRG(self, dmrg_input)
            chain_length = self.model.L;
            m_warmup = dmrg_input.m_warmup;
            conserved_QNum = self.model.conserved_QNum;
            
            assert(self.model.L >= 5, "Chain length must be greater than 4");
            
            %iDMRG step and save each system block corresponding to a
            %certain chain length to left block list and right block list
            fprintf("------------------------------------------------- \n");
            fprintf("Performing build up using iDRMG: \n");
            sys_block = self.model.InitSingleSiteBlock();         
            iteration_count = floor((chain_length-2)/2);
            direction = 'r';
            for iteration=1: iteration_count 
                env_block = sys_block;
                self.l_block_list(sys_block.length) = sys_block;
                self.r_block_list(env_block.length) = env_block;
                [sys_block, gnd_state, energy, rbi, T_matrix] = ... 
                    self.SingleDMRGStep(sys_block, env_block, m_warmup, direction, conserved_QNum);
            end 
            fprintf("infinite DRMG success with L=%d, E=%d \n", chain_length, energy);
            fprintf("\nPerforming fDMRG sweep: \n");
            env_block = sys_block;
            
            self.l_block_list(sys_block.length) = sys_block;
            self.l_T_matrices(sys_block.length) = T_matrix;      
            self.r_block_list(env_block.length) = env_block;   
            self.r_T_matrices(env_block.length) = T_matrix;
            
            if mod(chain_length,2) ~= 0
                [env_block, gnd_state, energy, rbi, T_matrix] = ...
                    self.SingleDMRGStep(sys_block, env_block, m_warmup, direction, conserved_QNum);  
                self.r_block_list(env_block.length) = env_block;
                self.r_T_matrices(env_block.length) = T_matrix;  
            end
            
            %sweeping for the finite DMRG 
            if isfield(dmrg_input, "sweep_list")
                sweep_list = dmrg_input.sweep_list;
                for i=1: length(sweep_list)
                    sweep_list = dmrg_input.sweep_list;
                    [sys_block,energy, gnd_state, rbi] = self.SingleSweep(sys_block, sweep_list(i), i);
                end
            elseif isfield(dmrg_input, "tolerence")
                m = m_warmup;
                tolerence = dmrg_input.tolerence;
                prev_energy = energy;
                [sys_block, energy, gnd_state, rbi] = self.SingleSweep(sys_block, m, 1);
                
                sweep = 1;
                while abs(energy-prev_energy) > tolerence
                    prev_energy = energy;
                    m = m + 5;
                    sweep = sweep + 1;
                    [sys_block,energy, gnd_state, rbi] = self.SingleSweep(sys_block, m, sweep);
                end
            end
            
            fprintf("Finite DRMG success with L=%d, E=%d \n", chain_length, energy);
            fprintf("------------------------------------------------- \n");
        end
%--------------------------------------------------------------------------        
        function [sys_block, energy, gnd_state, rbi] = SingleSweep(self, sys_block, m, sweep)
            conserved_QNum = self.model.conserved_QNum;
            sys_orig_length = sys_block.length;
            L = self.model.L;
            sys_label = 'l';
            env_label = 'r';
            while true
                if env_label == 'r'
                    env_block = self.r_block_list(L-sys_block.length-2); 
                else
                    env_block = self.l_block_list(L-sys_block.length-2); 
                end

                if env_block.length == 1
                    [sys_block, env_block] = deal(env_block, sys_block);
                    [sys_label, env_label] = deal(env_label, sys_label);    
                end

                direction = env_label;
                [sys_block, gnd_state, energy, rbi, T_matrix] = ... 
                    self.SingleDMRGStep(sys_block, env_block, m, direction, conserved_QNum);

                if sys_label == 'r'
                    self.r_block_list(sys_block.length) = sys_block;
                    self.r_T_matrices(sys_block.length) = T_matrix;
                else
                    self.l_block_list(sys_block.length) = sys_block;
                    self.l_T_matrices(sys_block.length) = T_matrix;
                end
                if (sys_label == 'l') && sys_block.length == sys_orig_length
                    fprintf("Energy found for sweep %d with m=%d: %d\n", sweep, m, energy);
                    self.energy = energy;
                    self.gnd_state = gnd_state;
                    self.rbi = rbi;
                    break;
                end
            end            
        end     
%--------------------------------------------------------------------------        
        function returned_measurements = Measurements(self, measurements)
            measurementContainer = MeasurementContainer(self, measurements);
            returned_measurements = measurementContainer.returned_measurements;
        end
    end
end

function map = KeyValueMap(array)
    map = containers.Map('KeyType','double','ValueType','any');
    possible_keys = unique(array);
    for i=1: length(possible_keys)
        possible_key = possible_keys(i);
        map(possible_key) = find(array == possible_key);
    end
end