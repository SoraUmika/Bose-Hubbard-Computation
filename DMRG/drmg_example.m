N = 3;
BoseHubbard = BoseHubbardChain(N,5,1,1,1);

model = BoseHubbard;
chain_length = 20;
m = 10;
target_QNum = N;
Gs_Energy = iDMRG(model, chain_length, m, target_QNum);

for M=4:2:chain_length
   H = BoseHubbard.Exact_H(M);
   Energy = eigs(H, 1, 'smallestreal');
   fprintf("L=%d, E=%d\n", M, Energy);
end