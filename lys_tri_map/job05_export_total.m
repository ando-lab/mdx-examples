%% Export diffuse map to h5 file

%% design grid
load('proc/unitCellInventory.mat','Crystal');
grid = load('proc/7_1/grid.mat','options');
ndiv = grid.options.ndiv;
smax = 0.8;

MT = proc.script.MapTools(...
    'SpaceGroup',symm.SpaceGroup(Crystal.spaceGroupNumber),...
    'Basis',latt.Basis(Crystal.a,Crystal.b,Crystal.c,Crystal.alpha,Crystal.beta,Crystal.gamma).invert,...
    'Grid',latt.PeriodicGrid(ndiv,[0,0,0],[1,1,1]),...
    'type','intensity');

%MT = MT.resize('radius',smax).resize('asu');
MT = MT.resize('roi',[-22.49,22.49,-26.49,26.49,-.49,28.49]);

%%
clear Crystal grid

%% map data to grid

load proc/mergeFine.mat
hklMerge.h = hklMerge.h + hklMerge.dh/ndiv(1);
hklMerge.k = hklMerge.k + hklMerge.dk/ndiv(2);
hklMerge.l = hklMerge.l + hklMerge.dl/ndiv(3);
hklMerge = hklMerge(:,{'h','k','l','I','sigma'});

hklMerge = hklMerge(~isinf(hklMerge.sigma),:); 

[I,sigma] = MT.table2array(hklMerge,'direct','replace',NaN);

msk = MT.spherical_mask(smax);
I(~msk) = NaN;
sigma(~msk) = NaN;
%%
MT.export('export/lys_tri_map.h5','total',I,sigma,'Datatype','single','ChunkSize',ndiv);
%%
load proc/mergeFineSplit.mat

hklMerge.h = hklMerge.h + hklMerge.dh/ndiv(1);
hklMerge.k = hklMerge.k + hklMerge.dk/ndiv(2);
hklMerge.l = hklMerge.l + hklMerge.dl/ndiv(3);
hklMerge = hklMerge(:,{'h','k','l','I1','sigma1','I2','sigma2'});

isIncl = ~isinf(hklMerge.sigma1) & ~isinf(hklMerge.sigma2);
hklMerge = hklMerge(isIncl,:); 

[I1,sigma1,I2,sigma2] = MT.table2array(hklMerge,'direct','replace',NaN);

I1(~msk) = NaN;
sigma1(~msk) = NaN;
I2(~msk) = NaN;
sigma2(~msk) = NaN;

MT.export('export/lys_tri_map.h5','total',I1,sigma1,I2,sigma2,'append',true,'Datatype','single','ChunkSize',ndiv);
