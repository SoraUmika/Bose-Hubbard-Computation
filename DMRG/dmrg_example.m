N = 60;
mu = 1;
t = 1;
U = 1;
BoseHubbard = BoseHubbardChain(N,U,mu,t,BoundCond.open);

model = BoseHubbard;
chain_length = 8;
m_warmup = 20;
m = 40;
target_QNum = N;
sweep_count = 3;

tic
BoseHubbard_DMRG = DMRG(model);
BoseHubbard_DMRG.fDMRG(chain_length, sweep_count, m_warmup, m, target_QNum);
toc

%test_exact_diagonalization(model, chain_length)
function none = test_exact_diagonalization(model, chain_length)
    fprintf("Exact Diagonalization: \n")
    M = chain_length;
    for i=4:2:M
        basis = model.Generate_Basis(i);
        [Gs, Energy] = model.ExactGsEnergy(i);
        fprintf("L=%d, E=%d\n", i, Energy); 
    end
     sqr_amp = Gs.^2;
     [max_comp, index_max] = max(sqr_amp);
     [min_comp, index_min] = min(sqr_amp);

    norm_basis = basis; 
    for i=1: size(basis,1)
        norm_basis(i, :) = norm_basis(i, :)*sqr_amp(i);
    end

    atom_dist = sum(norm_basis, 1);

    most_probable_eState = basis(index_max,:);
    least_probable_eState = basis(index_min,:);
    none = 1;
    %bar(atom_dist);
end
