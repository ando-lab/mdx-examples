%% compute variational intensity and add to export

mapFileName = 'export/lys_tet_map.h5';
slim = 0.625;
edges = 0.035:.01:0.625;

%% interpolate half-integer map and calculate intensity stats

%%
% design the target grid
MT = proc.script.MapTools.import(mapFileName,'total');

MI = proc.script.MapInterp(...
    'MT',MT,...
    'target_supercell',[2,2,2],...
    'target_slim',slim);
[To,MTo] = MI.target_hkl();

% select the half-integer points only
is_half_integer = @(x) logical(mod(round(x*2),2));
To = To(is_half_integer(To.h) & is_half_integer(To.k) & is_half_integer(To.l),:);

%%
% interpolate I, sigma maps
[I,sigma] = MT.read_data(mapFileName,'total','I','sigma');
[To] = MI.sgolay_filter_at_target(I,sigma,To,1/150);

[I,sigma] = MT.read_data(mapFileName,'total','I1','sigma1');
[To1] = MI.sgolay_filter_at_target(I,sigma,To,1/150);

[I,sigma] = MT.read_data(mapFileName,'total','I2','sigma2');
[To2] = MI.sgolay_filter_at_target(I,sigma,To,1/150);

clear I sigma
%%
% calculate intensity statistics

[I,sigma] = MTo.table2array(To,'direct','replace',NaN);
[I1,sigma1] = MTo.table2array(To1,'direct','replace',NaN);
[I2,sigma2] = MTo.table2array(To2,'direct','replace',NaN);

I(isinf(sigma)) = NaN;
I1(isinf(sigma1)) = NaN;
I2(isinf(sigma2)) = NaN;

SVR = proc.script.StatisticsVsRadius(...
        'edges',edges,...
        'Basis',MTo.Basis,...
        'PeriodicGrid',MTo.Grid);

stats = SVR.run(I); % calculate mean and variance
stats12 = SVR.run(I1,I2); % calculate mean, variance, covariance, and cc12

%%
% define the isotropic intensity component

cc12 = stats12.cc;
ccstar = sqrt(2*cc12./(1+cc12));

Iiso = stats.av - sqrt(stats.var).*ccstar;

%figure;plot(stats.r,Iiso)

%%
% make a smooth version of the isotropic part

IL = proc.scale.InterpLin1('Nw',101,'xmin',edges(1),'xmax',edges(end),'x',stats.r);

% fit the weighted mean
lambda = 1E1; % regularization parameter smoothness
pu0 = (IL.A'*IL.A + lambda*(IL.B'* IL.B))\(IL.A'*Iiso);

%% compute variational intensity map


[I,sigma] = MT.read_data(mapFileName,'total','I','sigma');

%% 
% fill in missing data points using nearest neighbors

isOK = ~isnan(I);

ker = zeros(3,3,3);
ker([1,3],2,2) = 1;
ker(2,[1,3],2) = 1;
ker(2,2,[1,3]) = 1;

A0 = I;
A0(~isOK) = 0;
n0 = convn(isOK,ker,'same');
A0 = convn(A0,ker,'same')./n0;
isnew = ~isOK & n0>0;
I(isnew) = A0(isnew);
sigma(isnew) = Inf; % signals that the value was interpolated

%% 
% subtract background

% compute s
[x,y,z] = MT.Grid.grid();
[x,y,z] = MT.Basis.frac2lab(x,y,z);
s = sqrt(x.*x + y.*y + z.*z);

% included points
isIncl = s <= IL.xmax & s >= IL.xmin & MT.isASU();

IL.x = s(isIncl);
Ibkg = NaN*ones(MT.Grid.N);
Ibkg(isIncl) = IL.interp(pu0);

I = I - Ibkg;
sigma(~isIncl) = NaN;

ismissing = isnan(I) & isIncl;
I(ismissing) = 0;
sigma(ismissing) = Inf;
%% export
ndiv = MT.Grid.invert.P;
MT.export(mapFileName,'variational',I,sigma,'newmap',true,'Datatype','single','ChunkSize',ndiv);
