%% export simulated scattering to h5 file
slim = 0.625;

%% load goodvibes simulation
[MTi,I] = proc.script.MapTools.import('export/lys_tet_goodvibes.h5','supercell','I');
tmp = MTi.array2table(I);
MTi = MTi.resize('symexpand');
I = MTi.table2array(tmp,'symexpand','replace',NaN);

%% load experimental data
mapFileName = '../lys_tet_map/export/lys_tet_map.h5';
[MTe,Ie,sigmaIe] = proc.script.MapTools.import(mapFileName,'total','I','sigma');

isIncl = ~isnan(Ie) & ~isnan(sigmaIe) & ~isinf(sigmaIe);

T = MTe.array2table(Ie,isIncl);
tmp = MTe.array2table(sigmaIe,isIncl);
T.sigma = tmp.(4);

[Ie,sigmaIe] = MTi.table2array(T,'symexpand','replace',NaN);
sigmaIe(isnan(Ie)) = Inf;
clear MTe isIncl T tmp

%% design a target grid

MI = proc.script.MapInterp(...
    'MT',MTi,...
    'target_supercell',[1,1,2],...
    'target_slim',slim + 0.05);

[To,MTo,msk] = MI.target_hkl();

%% run the filter routine

% sgolay filter (first pass)
Tfilt = MI.sgolay_filter_at_target(Ie - I,sigmaIe,To);

% interpolate sgolay result onto fine grid
I_sub_filt_interp = MI.interpolate_from_table(Tfilt); 

% sgolay filter (second pass) with robust bisquare weights
resid = Ie - I - I_sub_filt_interp;
[Tfilt,wbsq] = MI.sgolay_filter_at_target(Ie - I,sigmaIe,To,[],[],resid);

% symmetry expand bisquare weights
wbsq = MTi.table2array(MTi.array2table(wbsq),'symexpand','mean');
wbsq(isnan(wbsq)) = 0;

% symmetry expand filtered intensity
[I_sub_filt,sigmaI_sub_filt] = MTo.table2array(Tfilt,'symexpand');
sigmaI_sub_filt(isnan(I_sub_filt)) = Inf;
I_sub_filt(isnan(I_sub_filt)) = 0;

% compute the bisquare-weighted average signal
I_sub_filt_interp = MI.interpolate_from_table(Tfilt); 

I_sub_fill = Ie - I;
I_sub_fill(isnan(I_sub_fill)) = 0;

I_sub_fill = (1-wbsq).*I_sub_filt_interp + wbsq.*I_sub_fill;
sigmaI_sub_fill = sigmaIe./wbsq;

%% export subtracted / filled map

[MTout,cropfun] = MTi.resize('asu');

I = cropfun(I_sub_fill);
sigma = cropfun(sigmaI_sub_fill);

msk = MTout.spherical_mask(slim) & MTout.isASU();
I(~msk) = NaN;
sigma(~msk) = NaN;

h5Out = 'export/lys_tet_map_sub.h5';
chunksize = MTi.Grid.invert.P; 

MTout.export(h5Out,'sub_fill',I,sigma,'Datatype','single','ChunkSize',chunksize);

%% export filtered maps

%msk = MTo.spherical_mask(slim);% & MTo.isASU();

I = I_sub_filt;
sigma = sigmaI_sub_filt;

%I_sub(~msk) = NaN;
%sigma(~msk) = NaN;

MTo.export(h5Out,'sub_filt',I,sigma,'Datatype','single','newmap',true);


