clear all;
%First we specify the model parameters: N=NUM ATOMS; L=CHAIN LENGTH, 
%mu=CHEMICAL POTENTIAL, t=HOPPING CONSTANT, U=ON-SITE INTERACTION CONSTANT,
%V=NEAREST NEIGHBOR INTERACTION CONSTANT, k=HARMONIC-POTENTIAL CONSTANT
N = 6; L = 6; mu = 0; t = 1; U = 1; V = 0; k = 0; conserved_QNum = N;

%Input these parameters into the BoseHubbardChain class, it would then
%return a object associated with the paraters, think of this as a struct.
%This will be our model for DMRG
model = BoseHubbardChain(L, N, U, mu, t, V, k, BoundCond.open, conserved_QNum);

%we can run exact diagonalization as follow where we reference the model we
%just created, this allows the fucntion to use whatever parameters that was
%in the "model" we created
%[Energy, gnd_state] = model.ExactGsEnergy(); 

%% 
%First, create an DMRG object, this object takes in the model that we just
%created and can now be use to perform DMRG
BoseHubbard_DMRG = DMRG(model);

%There are two type of DMRG we can perform: iDMRG and fDMRG. Both of these
%takes in an structure named "dmrg_input". You need to pass in different
%field parameters depending on which method you use. For iDMRG, you only
%need "m_warmup" parameter, this will be the maximum number of basis that
%you kept in the DMRG step. Since fDMRG builds on top of iDRMG, it is
%necessary to include m_warmup in fDMRG as well.
%dmrg_input.m_warmup = 80;
%[Energy2, gnd_state] = BoseHubbard_DMRG.iDMRG(dmrg_input);

%For fDMRG, we sweep through the chain and gains better accuracy by using
%better and better values of "m". Here we provide two options for this. The
%first option is manually adding a sweep list that contains increasing values
%of "m" that you wish to perform the sweep. You don't want to be using the
%same "m" value for every sweep because generally your accuracy no longer 
%increases after just a single sweep using a particular value of "m"

dmrg_input.m_warmup = 10;
dmrg_input.sweep_list = [10,20, 30, 40];
[Energy2, gnd_state] = BoseHubbard_DMRG.fDMRG(dmrg_input);

%The second option is to provide a "tolerence" field that will determine
%the accuracy of the energy that you wish to have. For this method, we start
%out with some default value of "m" and increasing it by some value for each 
%additional sweep. For each additional sweep, we compare the
%difference in energies between each successive sweep with the tolerence.
%The smaller the tolernce, the more time it will take and the more accurate
%it will be. However, although this method is convenient, there is some
%subtleties to it.
%dmrg_input.m_warmup = 10;
%dmrg_input.tolerence = 1e-3;
%[Energy, gnd_state] = BoseHubbard_DMRG.fDMRG(dmrg_input);


%To perform measurements, you must use fDMRG, in addition you must specify 
%which measurements you would like to take. Following is an example of
%measurements of n_i and <n_i n_j>. In general, if you want to perform an
%measurement on a particular site, then in your model, in this case the
%bosehubabrd model must include the single site operator corresponding to
%your measurement. For instance, what is already included in our single
%site operator is "n" and "b" and hence it is able to perform on site
%measurements for these operators. 
measurements = struct('site_indices',{},'operator_names',{});
for i=1: L
    measurement.operator_names = 'n';
    measurement.site_indices = i;
    measurements(i) = measurement;
end
for i=1: L
    for j=1:L
        measurement.operator_names = ['n','n'];
        measurement.site_indices = [i,j];
        measurements(end+1) = measurement;
    end
end

returned_measuremetns = BoseHubbard_DMRG.Measurements(measurements);


