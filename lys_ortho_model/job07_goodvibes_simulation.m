%% simulate 1-phonon diffuse scattering from refined lattice dynamics model

load proc/atomic_model.mat Basis SpaceGroup
[MT,rho] = proc.script.MapTools.import('export/lys_ortho_edens.h5','unwrapped_cropped','rho');

load proc/goodvibes_model.mat ENM
load proc/goodvibes_fit.mat LD

supercell = LD.supercell;
V = LD.V;

LDT = proc.script.LatticeDynamicsTools('supercell',supercell,'Cell',ENM.Cell,'M',LD.M);
%LDT.interp_osr = 3; % to avoid memory issues on home computer.
clear ENM LD

%% simulation grid (use experimental map)
mapFileName = '../lys_ortho_map/export/lys_ortho_map.h5';
slim = 0.625;
MTsim = proc.script.MapTools.import(mapFileName,'total');
isIncl = MTsim.isASU() & MTsim.spherical_mask(slim);

[h,k,l] = MTsim.Grid.grid();
h = h(isIncl); k = k(isIncl); l = l(isIncl);

%% compute the one phonon structure factors by interpolation

sffun = LDT.calc1PSFInterpFromMap(rho,MT.Grid,MT.Basis);

%% simulate

I = NaN*ones(MTsim.Grid.N);
I(isIncl) = LDT.calc1PIntensity(V,sffun,h,k,l);

%% export the results

chunksize = MTsim.Grid.invert.P; % = 1/P.delta = [5,5,11]

MTsim.export('export/lys_ortho_goodvibes.h5','supercell',I,'Datatype','single','ChunkSize',chunksize);
