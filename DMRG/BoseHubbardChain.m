classdef BoseHubbardChain
    properties
        t
        U
        mu
        N
        bound_cond
        d
        single_site_QNum
        n_op
        b_op
        H_op 
    end
    
    methods
        function self = BoseHubbardChain(N, U, mu, t, bound_cond)           
            %Bose Hubbard Model Parrameters
            self.t = t;
            self.U = U;
            self.mu = mu;
            
            %open_bound=1 => open boundary, open_bound=0 => closed boundary
            self.bound_cond = bound_cond;
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
            self.H_op = diag(U/2 * ndiag .* (ndiag - 1) - mu * ndiag);            
        end
        
        %A block consists of four parameters: operator list, basis_size, 
        %length, and possible quantum numbers with each index "i" marking the
        %ith basis vector. 
        
        %here we initialize the block for a single site
        function block = InitSingleSiteBlock(self)
            block.op_list = containers.Map;
            if self.bound_cond == BoundCond.open
                block.op_list('H') = sparse(self.H_op); 
                block.op_list('conn_n') = sparse(self.n_op); 
                block.op_list('conn_b') = sparse(self.b_op); 
            elseif self.bound_cond == BoundCond.periodic
                block.op_list('H') = sparse(self.H_op); 
                block.op_list('l_n') = sparse(self.n_op);
                block.op_list('r_n') = sparse(self.n_op);
                block.op_list('l_b') = sparse(self.b_op);
                block.op_list('r_b') = sparse(self.b_op);
            else
                fprintf("Error: No valid boundary condition! \n");
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
        function enlarged_block = EnlargeBlock(self, block, direction)
            new_op_list = containers.Map;
            block_d = block.basis_size;
            if self.bound_cond == BoundCond.open
                new_op_list('H') = sparse(kron(block.op_list('H'), eye(self.d)) ...
                    + kron(eye(block_d), self.H_op) ...
                    + self.SiteSite_H(block.op_list('conn_b'),self.b_op));
                new_op_list('conn_b') = sparse(kron(eye(block_d), self.b_op));
                new_op_list('conn_n') = sparse(kron(eye(block_d), self.n_op));
            elseif self.bound_cond == BoundCond.periodic
                assert(direction == 'l' | direction == 'r');
                if direction == 'l'
                    new_op_list('H') = sparse(kron(block.op_list('H'), eye(self.d)) ...
                        + kron(eye(block_d), self.H_op) ...
                        + self.SiteSite_H(block.op_list('l_b'), self.b_op));     
                    new_op_list('l_n') = sparse(kron(eye(block_d), self.n_op));
                    new_op_list('l_b') = sparse(kron(eye(block_d), self.b_op));
                    new_op_list('r_n') = sparse(kron(block.op_list('r_n'), eye(self.d)));
                    new_op_list('r_b') = sparse(kron(block.op_list('r_b'), eye(self.d)));
                else
                    new_op_list('H') = sparse(kron(block.op_list('H'), eye(self.d)) ...
                        + kron(eye(block_d), self.H_op) ...
                        + self.SiteSite_H(block.op_list('r_b'),self.b_op)); 
                    new_op_list('l_n') = sparse(kron(block.op_list('l_n'), eye(self.d)));
                    new_op_list('l_b') = sparse(kron(block.op_list('l_b'), eye(self.d)));   
                    new_op_list('r_n') = sparse(kron(eye(block_d), self.n_op));
                    new_op_list('r_b') = sparse(kron(eye(block_d), self.b_op));               
                end              
            end
            enlarged_block.op_list = new_op_list;
            enlarged_block.basis_size = block_d*self.d;
            enlarged_block.length = block.length + 1;
            basis_QNum = outerop(block.basis_QNum ,self.single_site_QNum, '+');
            enlarged_block.basis_QNum = reshape(basis_QNum',1,[]);
        end
        
        %construct super block given a system block (or left block) and
        %a environment block (or right block)
        function superBlock_H = ConstructSuperBlock_H(self, sys_block, env_block)
            if self.bound_cond == BoundCond.open
               superBlock_H = kron(sys_block.op_list('H'), eye(env_block.basis_size)) + ...
                   kron(eye(sys_block.basis_size), env_block.op_list('H')) + ...
                   self.SiteSite_H(sys_block.op_list('conn_b'), env_block.op_list('conn_b'));
            elseif self.bound_cond == BoundCond.periodic
               superBlock_H = kron(sys_block.op_list('H'), eye(env_block.basis_size)) + ...
                   kron(eye(sys_block.basis_size), env_block.op_list('H')) + ...
                   self.SiteSite_H(sys_block.op_list('r_b'), env_block.op_list('l_b')) + ...
                   self.SiteSite_H(sys_block.op_list('l_b'),env_block.op_list('r_b'));
            end
        end
        
        function basis_vectors = Generate_Basis(self, M)
            D = factorial(self.N+M-1)/(factorial(self.N)*factorial(M-1));
            D = round(D);
            basis_vectors = zeros(D,M);
            basis_vectors(1,1) = self.N;
            current_row = 1;
            while basis_vectors(current_row, M) ~= self.N
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
                basis_vectors(next_row, k+1) = self.N - sum(basis_vectors(next_row, 1:k));
                current_row = next_row; 
            end
        end 
        
        function [Gs, Energy] = ExactGsEnergy(self, M)
            basis = self.Generate_Basis(M);
            D = factorial(self.N+M-1)/(factorial(self.N)*factorial(M-1));
            D = round(D);
            H = zeros(D);
            T = zeros(1,D);
            for i=1:D
              tag = Tag(basis(i,:));
              T(1,i) = tag;
            end  
            %[T,indices] = sort(T);
            %basis = basis(indices,:);
            for v=1:D
                v_vect.v = basis(v, :);
                v_vect.c = 1;      
                int_term = 0;
                chem_term = 0;
                for i=1:M
                    site = i;
                    site_adj = i+1;
                    
                    vect3 = Num_Op(site, v_vect);
                    n = vect3.c;
                    int_term = int_term + n*n-n;
                    chem_term = chem_term + n;      
                    
                    if site == M
                        if self.bound_cond == BoundCond.open
                            break;
                        elseif self.bound_cond == BoundCond.periodic
                            site_adj = 1;                            
                        end
                    end

                    vect1 = Create_Op(site_adj,Annih_Op(site, v_vect));
                    vect2 = Create_Op(site,Annih_Op(site_adj ,v_vect));
                    
                    tag1 = Tag(vect1.v);
                    tag2 = Tag(vect2.v);
                    index1 = find(T==tag1);
                    index2 = find(T==tag2);
                    H(index1, v) = H(index1, v) - self.t *vect1.c;
                    H(index2, v) = H(index2, v) - self.t *vect2.c;

                end 
                H(v,v) = self.U/2*int_term - self.mu*chem_term;  
            end
            Hamiltonian = sparse(H);
            [Gs, Energy] = eigs(Hamiltonian, 1, 'smallestreal'); 
            
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

function tag = Tag(basis_vect)
    M = length(basis_vect);
    p_vect = 100*(1:M)+3;
    tag = sum(basis_vect.*sqrt(p_vect),'all'); 
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
