%N is the number of spins of the chain
function Hamiltonian = Heisenberg_Hamiltonian(N)
    function Basis = Heisenberg_Basis(N)
        Basis = dec2bin(0:2^N-1)' - '0';
        Basis = transpose(Basis);
    end
    
    Basis = Heisenberg_Basis(N);
    %N is the site count and D the the total number of basis
    [D,N] = size(Basis);
    % the vector is in the form of (constant, vector_elements)
    % index refers to the site number

    function vect = S_z(index, basis_vect)
        vect.v = basis_vect.v;
        if basis_vect.v(index) == 1
            vect.c = 1/2*basis_vect.c;
        else
            vect.c = -1/2*basis_vect.c;
        end 
    end

    function vect = S_p(index, basis_vect)
        vect.v = basis_vect.v;
        if basis_vect.v(index) == 1
            vect.c = 0;
        else
            vect.c = 1*basis_vect.c;
            vect.v(index) = 1;
        end 
    end

    function vect = S_m(index, basis_vect)
        vect.v = basis_vect.v;
        if basis_vect.v(index) == 1
            vect.c = 1*basis_vect.c;
            vect.v(index) = 0;
        else
            vect.c = 0;
        end 
    end

    function value = Bra_Ket(vect1,vect2)
        if isequal(vect1.v,vect2.v)
            value = vect2.c*vect1.c;
        else
            value = 0;
        end
    end

    J = -2;
    Hamiltonian = zeros(D,D);
    for i=1: D
        for j=1: D
            for k=1: N-1
                basis_a.v = Basis(i,:);
                basis_a.c = 1;
                basis_b.v = Basis(j,:);
                basis_b.c = 1;
                index = k;
            
                value_1 = Bra_Ket(basis_a, S_z(index, S_z(index+1, basis_b)));
                value_2 = 0.5 * Bra_Ket(basis_a, S_p(index, S_m(index+1, basis_b)));
                value_3 = 0.5 * Bra_Ket(basis_a, S_m(index, S_p(index+1, basis_b)));
                kth_contribution = value_1 + value_2 +value_3;
                Hamiltonian(i,j) = Hamiltonian(i,j) + kth_contribution;      
            end
            Hamiltonian(i,j) = -J*Hamiltonian(i,j);
        end     
    end 
end

