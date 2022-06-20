function [Site_Num, Energies] = DMRG_Heisenberg(m,target_site_num)
    %initialize operators for single site
    I_d = eye(2);
    S_z = [1/2 0;0 -1/2];
    S_p = [0 1; 0 0];
    S_m = [0 0; 1 0];

    %We begin with one site and call this site block A
    S_zA = S_z;
    S_pA = S_p;
    S_mA = S_m;
    I_dA = I_d;
    H_A = zeros(2);
    J = -2;
    Energy = -0.75;
    
    %number of sites = 2*cycle_count+2
    cycle_count = (target_site_num-2)/2;
    Energies = zeros(cycle_count,1);
    Site_Num = zeros(cycle_count,1);
    for cycle=1: cycle_count
        %Now we add a site to block A, and we find the operators of block A + site
        H_A = kron(H_A,I_d) -J*(kron(S_zA, S_z) + 0.5*(kron(S_pA, S_m) + ...
            kron(S_mA, S_p))); 

        S_zA = kron(I_dA, S_z);
        S_pA = kron(I_dA, S_p);
        S_mA = kron(I_dA, S_m);
        I_dA = kron(I_dA, I_d);

        %Formulate the superblock using reflection symmetry
        H_superBlock = kron(H_A, I_dA) + kron(I_dA,H_A) - ...
            J*(0.5*(kron(S_pA, S_mA) + kron(S_mA, S_pA)) + kron(S_zA, S_zA));
        H_superBlock = 0.5*(H_superBlock+H_superBlock');

        %Getting the ground state from exat diagonalization
        %opts.disp = 0;  
        %opts.issym = 1;
        %opts.real = 1;
        LastEnergy = Energy;
        [Psi, Energy] = eigs(H_superBlock,1,'smallestreal');
        EnergyPerBond = (Energy - LastEnergy) / 2;
        
        Site_Num(cycle) = cycle*2+2;
        Energies(cycle) = EnergyPerBond;
        
        %Getting the reduce density matrix
        nr=size(Psi,1);
        Dim = sqrt(nr);
        PsiMatrix = reshape(Psi, Dim, Dim);
        Rho = PsiMatrix * PsiMatrix';

        %Diagonlize the reduce density matrix
        [V, D] = eig(Rho);
        [D, Index] = sort(diag(D), 'descend');  % sort eigenvalues descending
        V = V(:,Index);                     % sort eigenvectors the same way

        %Truncation of the number of basis vectors
        NKeep = min(size(D, 1), m); % numbers of basis vectors to keep
        T = V(:, 1:NKeep); % basis vectors kepted
        TruncationError = 1 - sum(D(1:NKeep)); 

        S_zA = T'*S_zA*T;
        S_pA = T'*S_pA*T;
        S_mA = T'*S_mA*T;
        I_dA = T'*I_dA*T;
        H_A = T'*H_A*T;
    end
end


