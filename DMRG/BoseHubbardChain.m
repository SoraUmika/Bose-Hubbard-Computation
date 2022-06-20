classdef BoseHubbardChain
    properties
        t
        U
        mu
        N
        open_bound
        d
        single_site_QNum
        n_op
        b_op
        H_op 
    end
    
    methods
        function self = BoseHubbardChain(N, U, mu, t, open_boundary_cond)           
            %Bose Hubbard Model Parrameters
            self.t = t;
            self.U = U;
            self.mu = mu;
            
            %open_bound=1 => open boundary, open_bound=0 => closed boundary
            self.open_bound = open_boundary_cond;
            self.N = N;
            
            %single site dimension
            self.d = N+1;
            ndiag = 0:self.d-1;
                            
            %Here nidag is an array, representing the possible particle 
            %number in a single site. For N particles, we have ndiag 
            %ranging from 0 to N. The location, or index in which these 
            %possible particle numbers are reprsents the basis number. 
            %So for example at index 5 if we have possible particle # is 4,
            %then we mean the 5th basis vector has a occupational number of
            %4 for the single site. 
            self.single_site_QNum = ndiag;
            
            %single site operators
            self.n_op = diag(ndiag);
            self.b_op = diag(sqrt(ndiag(2:self.d)),1);
            self.H_op = diag(0.5 * U * ndiag .* (ndiag - 1) - mu * ndiag);            
        end
        
        %A block consists of four parameters: operator list basis_size, 
        %length, and possible quantum numbers with each index "i" marking the
        %ith basis vector. 
        
        %here we initialize the block for a single site
        function block = InitSingleSiteBlock(self)
            block.op_list = containers.Map;
            if self.open_bound    
                block.op_list('H') = sparse(self.H_op); 
                block.op_list('conn_n') = sparse(self.n_op); 
                block.op_list('conn_b') = sparse(self.b_op); 
            else        
                block.op_list('H') = sparse(self.H_op); 
                block.op_list('l_n') = sparse(self.n_op); 
                block.op_list('l_b') = sparse(self.b_op); 
                block.op_list('r_n') = sparse(self.n_op); 
                block.op_list('r_b') = sparse(self.b_op); 
            end        
            block.basis_size = self.d;
            block.basis_QNum = self.single_site_QNum;
            block.length = 1;
        end
        
        %this calculates the hamiltonian of two adjacent sites
        function H = SiteSite_H(self,b_i,b_j)
            H = -self.t * (kron(b_i, b_j') + kron(b_i', b_j));
        end
        
        %enlarge the block by adding a single site by calculating the 
        %corresponding block parameters. 
        function enlarged_block = EnlargeBlock(self, block)
            enlarged_block.op_list = containers.Map;
            if self.open_bound
                enlarged_block.op_list('H') = sparse(kron(block.op_list('H'), eye(self.d)) ...
                    + kron(eye(block.basis_size), self.H_op) ...
                    + self.SiteSite_H(block.op_list('conn_b'),self.b_op));
                enlarged_block.op_list('conn_b') = sparse(kron(eye(block.basis_size), self.b_op));
                enlarged_block.op_list('conn_n') = sparse(kron(eye(block.basis_size), self.n_op));
            else 
                fprintf("closed boundary condition currently not avaiable");
                
            end
            enlarged_block.basis_size = block.basis_size*self.d;
            enlarged_block.length = block.length + 1;
            basis_QNum = outerop(block.basis_QNum ,self.single_site_QNum, '+');
            enlarged_block.basis_QNum = reshape(basis_QNum',1,[]);
        end
        
        %construct super block given a system block (or left block) and
        %a environment block (or right block)
        function superBlock_H = ConstructSuperBlock_H(self, sys_block, env_block)
            if self.open_bound
               superBlock_H = sparse(kron(sys_block.op_list('H'), eye(env_block.basis_size)) + ...
                   kron(eye(sys_block.basis_size), env_block.op_list('H')) + ...
                   self.SiteSite_H(sys_block.op_list('conn_b'), env_block.op_list('conn_b')));
            else
                fprintf("closed boundary condition currently not avaiable");
            end
        end
        
        %This function is currently not complete.
        function Energy = ExactGsEnergy(self, chain_length)
            currentBlock = self.InitSingleSiteBlock();
            Energy = eigs(currentBlock.op_list('H'),1,'smallestreal');
            
            for M=1: chain_length-1            
                D = factorial(self.N+M)/(factorial(self.N)*factorial(M));                     
                possible_QNum = unique(currentBlock.basis_QNum);
                restricted_indices = zeros(1,D);  
                new_basis_QNum = zeros(1,D);
                current_r_index = 1;
                for i=1: length(possible_QNum)
                    QNum = possible_QNum(i);
                    site_QNum = self.N - QNum;
                    site_QNum_index = find(self.single_site_QNum == site_QNum, 1);
                    if ~isempty(site_QNum_index) 
                        QNum_indices = find(currentBlock.basis_QNum == QNum);
                        for j=1:length(QNum_indices)
                          offset = currentBlock.basis_size*(QNum_indices(j)-1);
                          restricted_indices(current_r_index) = offset ...
                            + site_QNum_index;
                          new_basis_QNum(current_r_index) = QNum;
                          current_r_index = current_r_index + 1;                         
                        end
                    end                   
                end
                
                currentBlock = self.EnlargeBlock(currentBlock);
                
                currentBlock.basis_size = D;
                currentBlock.basis_QNum = new_basis_QNum;
                
                op_names = keys(currentBlock.op_list);              
                for i=1: currentBlock.op_list.Count
                    op_name = op_names(i);
                    op_matrix = currentBlock.op_list(op_name{1});

                    currentBlock.op_list(op_name{1}) = ...
                        op_matrix(restricted_indices,restricted_indices);

                end
                    
                Energy = eigs(currentBlock.op_list('H'),1,'smallestreal');              
                
                fprintf('L=%d, E=%d\n', M+1, Energy);                
            end
        end
        
        function Hamiltonian = Exact_H(self, M)
            D = factorial(self.N+M-1)/(factorial(self.N)*factorial(M-1));
            basis = Generate_Basis(self.N,M);
            Hamiltonian = zeros(D,D);
            for u=1:D
                for v=1:D
                    u_vect.c = 1;
                    u_vect.v = basis(u,:);
                    v_vect.c = 1;
                    v_vect.v = basis(v,:);       
                    element_kinetic = 0;
                    for i=1:M-1
                        vect1 = Create_Op(i+1,Annih_Op(i,v_vect));
                        vect2 = Create_Op(i,Annih_Op(i+1,v_vect));
                        element_kinetic = element_kinetic + Bra_Ket(u_vect,vect1) + Bra_Ket(u_vect, vect2);
                    end         
                    element_int = 0;
                    element_chem = 0;
                    for i=1:M
                        vect1 = Num_Op(i, Num_Op(i,v_vect));
                        vect2 = Num_Op(i, v_vect);
                        element_int = element_int + Bra_Ket(u_vect, vect1) - Bra_Ket(u_vect,vect2);
                        element_chem = element_chem + Bra_Ket(u_vect,vect2);
                    end
                    element = -self.t*element_kinetic + self.U/2*element_int ...
                        -self.mu*element_chem;       
                    Hamiltonian(u,v) = element;
                end
            end   
            Hamiltonian = sparse(Hamiltonian);
        end
    end
end

%imagine this function as an outer product, but instead of the operation of
%multiplication, we can also use the operation of addition.
function y=outerop(a,b,operator)
    if nargin<3
        operator='+';                     
    end  
    if isequal(operator,'*')               
        y=a(:)*b(:)';
    else    
      outera=a(:)*ones(1,length(b));       
      outerb=ones(length(a),1)*b(:).';     
      functionHandle=str2func(operator); 
      y=functionHandle(outera,outerb);    
    end
end


function basis_vectors = Generate_Basis(N,M)
    D = factorial(N+M-1)/(factorial(N)*factorial(M-1));
    basis_vectors = zeros(D,M);
    basis_vectors(1,1) = N;
    current_row = 1;
    while basis_vectors(current_row, M) ~= N
        k = 0;
        for j = M-1: -1: 1
            if basis_vectors(current_row, j) > 0
                k = j;
                break;
            end 
        end 
        next_row = current_row + 1;
        for j = 1: M
            if 1 <= j && j <= k-1
                basis_vectors(next_row,j) = basis_vectors(current_row,j);
            elseif j >= k+2
                basis_vectors(next_row,j) = 0;
            end
        end
        basis_vectors(next_row, k) = basis_vectors(current_row, k) - 1;
        basis_vectors(next_row, k+1) = N - sum(basis_vectors(next_row, 1:k));
        current_row = next_row; 
    end
end 

function result = Create_Op(index,vect)
    vect.v(index) = vect.v(index) + 1;
    vect.c = vect.c*sqrt(vect.v(index));
    result = vect;
end

function result = Annih_Op(index,vect)
    if vect.v(index) == 0
        vect.c = 0;
    else 
        vect.v(index) = vect.v(index)-1;
        vect.c = vect.c*sqrt(vect.v(index)+1);
    end
    result = vect;
end

function result = Num_Op(index, vect)
    vect.c = vect.c*vect.v(index);
    result = vect;
end  

function result = Bra_Ket(bra,ket)
    if ket.c == 0
        result = 0;
    elseif isequal(bra.v, ket.v)
        result = ket.c;
    else
        result = 0;
    end
end
