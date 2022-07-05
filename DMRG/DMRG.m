classdef DMRG
    properties
        model
        l_block_list 
        r_block_list 
    end
    
    methods
        function self = DMRG(model)
            self.model = model;
            self.l_block_list = containers.Map('KeyType','uint64','ValueType','any');
            self.r_block_list = containers.Map('KeyType','uint64','ValueType','any');
        end
        
        function [sys_block, gnd_state, energy] = ...
                SingleDMRGStep(self, sys_block, env_block, m, direction, conserved_QNum)                      
            sys_block = self.model.EnlargeBlock(sys_block, direction);
            env_block = self.model.EnlargeBlock(env_block, direction);                        
            superBlock_H = self.model.ConstructSuperBlock_H(sys_block, env_block);
            %basis_QNum is an arary that contains the value of all possible
            %quantum number for that block. The following KeyValueMap maps
            %the quantum numbers to the indices that they were found. For
            %example, consider the sys_QNum_Indices_Map. If the quantum number 
            %3 is found at indices 3,7,2 then sys_QNum_Indices_Map(3)
            %will give us {[3,7,2]}, note that this is in a cell array. 
            sys_QNum_Indices_Map = KeyValueMap(sys_block.basis_QNum);
            env_QNum_Indices_Map = KeyValueMap(env_block.basis_QNum);
            
            restricted_basis_indices = [];       
            
            %Now we consider the restriction of conserved quantum numbers
            sys_possible_QNum_List = cell2mat(keys(sys_QNum_Indices_Map));
            sys_QNum_Count = sys_QNum_Indices_Map.Count;
            new_sys_QNum_Indices_Map = ... 
                containers.Map(sys_possible_QNum_List, cell(1,sys_QNum_Count)); 
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
                            restricted_basis_indices(end+1) = index_offset + env_QNum_Indices(k);                   
                            new_index = length(restricted_basis_indices);                      
                            tmp_arr = new_sys_QNum_Indices_Map(sys_QNum);
                            new_sys_QNum_Indices_Map(sys_QNum) = [tmp_arr new_index];
                        end
                    end
                end
            end
            restricted_superBlock_H = ... 
                superBlock_H(restricted_basis_indices,restricted_basis_indices);
            [restricted_psi0, energy] = eigs(restricted_superBlock_H,1,'smallestreal');
            gnd_state = restricted_psi0;

            new_sys_QNum = keys(new_sys_QNum_Indices_Map);
            new_sys_Indices = values(new_sys_QNum_Indices_Map);   

            %we now calculate the reduce density matrix. Here, we make them in
            %block matrix form where each block of the matrix represents the 
            %corresponding quantum numbers 
            %rho_e_values = zeros(1,sys_block.basis_size);
            rho_e_values = [];
            rho_e_vects = zeros(sys_block.basis_size,1);   
            new_basis_QNum = [];
            current_evalue_index = 1;
            for i=1:new_sys_QNum_Indices_Map.Count
                if ~isempty(new_sys_Indices{i})
                    psi0_QNum = restricted_psi0(new_sys_Indices{i});          
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
            %reorder the eigralues, eigen vectors, the corresponding quantum 
            %number in descending eigenvalue order. 
            [rho_e_values, Index] = sort(rho_e_values, 'descend');
            rho_e_vects = rho_e_vects(:,Index);        
            new_basis_QNum = new_basis_QNum(Index);

            %truncation 
            NKeep = min(length(rho_e_values), m);
            T = rho_e_vects(:, 1:NKeep);
            new_basis_QNum = new_basis_QNum(:, 1:NKeep);
            trunc_err = 1 - sum(rho_e_values(1:NKeep)); 
            
            total_L = sys_block.length+env_block.length;
            fprintf("L=%i, E=%d, err=%d\n", total_L, energy, trunc_err);
            
            %update our curent block
            newBlock.basis_size = NKeep;
            newBlock.length = sys_block.length;
            newBlock.op_list = containers.Map;
            newBlock.basis_QNum = new_basis_QNum;
            sys_op_name = keys(sys_block.op_list); 
            for i=1: length(sys_op_name)
               op_name = sys_op_name(i);
               op_mat = sys_block.op_list(op_name{1});
               newBlock.op_list(op_name{1}) = T'*op_mat*T;
            end
            
            sys_block = newBlock;                  
        end
        
        function [gnd_state, energy] = iDMRG(self, chain_length, m, conserved_QNum)
            l_block = self.model.InitSingleSiteBlock();   
            iteration_count = (chain_length-2)/2;
            direction = 'r';
            for iteration=1: iteration_count 
                r_block = l_block;
                [l_block, gnd_state, energy] = ...
               self.SingleDMRGStep(l_block, r_block, m, direction, conserved_QNum);
            end
        end
        
        function [gnd_state, energy] = fDMRG(self, chain_length, sweep_count, m_warmup, m, conserved_QNum)
            fprintf("Performing build up using iDRMG: \n");
            sys_block = self.model.InitSingleSiteBlock();         
            iteration_count = (chain_length-2)/2;
            for iteration=1: iteration_count 
                env_block = sys_block;
                self.l_block_list(sys_block.length) = sys_block;
                self.r_block_list(env_block.length) = env_block;
                [sys_block, gnd_state, energy] = ... 
                    self.SingleDMRGStep(sys_block, env_block, m_warmup, 'r', conserved_QNum);
            end
            
            fprintf("\nPerforming fDMRG sweep: \n");
            env_block = sys_block;
            self.l_block_list(sys_block.length) = sys_block;
            self.r_block_list(env_block.length) = env_block;            
            
            for sweep=1: sweep_count
                fprintf("Sweep #: %d\n", sweep);
                L = chain_length;
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
                    [sys_block, gnd_state, energy] = ... 
                        self.SingleDMRGStep(sys_block, env_block, m, direction, conserved_QNum);

                    if sys_label == 'r'
                        self.r_block_list(sys_block.length) = sys_block;
                    else
                        self.l_block_list(sys_block.length) = sys_block;
                    end

                    if sys_label == 'l' && (2 * sys_block.length) == L
                        break;
                    end
                end            
            end
            
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