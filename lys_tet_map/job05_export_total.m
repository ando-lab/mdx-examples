%% Export diffuse map to h5 file

%% design grid
load('proc/unitCellInventory.mat','Crystal');
grid = load('proc/2/grid.mat','options');
ndiv = grid.options.ndiv;
smax = 0.625;
mapFileName = 'export/lys_tet_map.h5';

MT = proc.script.MapTools(...
    'SpaceGroup',symm.SpaceGroup(Crystal.spaceGroupNumber),...
    'Basis',latt.Basis(Crystal.a,Crystal.b,Crystal.c,Crystal.alpha,Crystal.beta,Crystal.gamma).invert,...
    'Grid',latt.PeriodicGrid(ndiv,[0,0,0],[1,1,1]),...
    'type','intensity');

MT = MT.resize('roi',[-.49,50.49,-.49,50.49,-.49,24.49]);

clear Crystal grid

%% map data to grid

load proc/mergeFine.mat
hklMerge.h = hklMerge.h + hklMerge.dh/ndiv(1);
hklMerge.k = hklMerge.k + hklMerge.dk/ndiv(2);
hklMerge.l = hklMerge.l + hklMerge.dl/ndiv(3);
hklMerge = hklMerge(:,{'h','k','l','I','sigma'});

hklMerge = hklMerge(~isinf(hklMerge.sigma),:); 

[I,sigma] = MT.table2array(hklMerge,'symexpand','replace',NaN);

msk = MT.spherical_mask(smax) & MT.isASU();
I(~msk) = NaN;
sigma(~msk) = NaN;
%%
MT.export(mapFileName,'total',I,sigma,'Datatype','single','ChunkSize',ndiv);
%%
load proc/mergeFineSplit.mat

hklMerge.h = hklMerge.h + hklMerge.dh/ndiv(1);
hklMerge.k = hklMerge.k + hklMerge.dk/ndiv(2);
hklMerge.l = hklMerge.l + hklMerge.dl/ndiv(3);
hklMerge = hklMerge(:,{'h','k','l','I1','sigma1','I2','sigma2'});

isIncl = ~isinf(hklMerge.sigma1) & ~isinf(hklMerge.sigma2);
hklMerge = hklMerge(isIncl,:); 

[I1,sigma1,I2,sigma2] = MT.table2array(hklMerge,'symexpand','replace',NaN);

I1(~msk) = NaN;
sigma1(~msk) = NaN;
I2(~msk) = NaN;
sigma2(~msk) = NaN;

MT.export(mapFileName,'total',I1,sigma1,I2,sigma2,'append',true,'Datatype','single','ChunkSize',ndiv);
