%% Lattice Dynamics Model for Tetragonal Lysozyme


load proc/atomic_model.mat Atoms SpaceGroup Basis

ENT = proc.script.ElasticNetworkTools.initialize(Atoms,Basis,SpaceGroup);
[g,Tg] = ENT.find_atom_groups();
ENT.groupAttributes = {}; % don't split into subgroups

[T,G,index] = ENT.externalModel();
ENM = ENT.exportModel(T,G);
M = ENT.M;

ResidueGroups = Tg(g(index{1}),:);

%% save result
save proc/goodvibes_model.mat ENM M ResidueGroups

