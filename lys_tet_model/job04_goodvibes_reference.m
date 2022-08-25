%% fetch reference halos and pre-compute structure factor amplitudes

%% reference halos to fit

% get the 400 most-intense Bragg reflections in this resolution bin
numRef = 400;

% choose reflections between 2.5 and 2.0A resolution
resMin = 2;
resMax = 2.5;

load proc/atomic_model.mat hklTable Basis

hklRef = selectStrongBraggPeaks(hklTable,Basis,numRef,resMin,resMax);
%% 
% load diffuse data around those peaks

mapFileName = '../lys_tet_map/export/lys_tet_map.h5';
dset = '/maps/variational';

[hklGrid,I,sigma] = proc.script.LatticeDynamicsTools.loadBrillouinZones(mapFileName, dset, hklRef);

%%
% save the reference data

save('proc/reference_halos.mat','hklGrid','I','sigma');

%% one-phonon structure factor calculations for reference data points

%% 
% calculate one-phonon structure factors from electron density

[MT,rho] = proc.script.MapTools.import('export/lys_tet_edens.h5','unwrapped_cropped','rho');
load proc/goodvibes_model.mat ENM ResidueGroups M

%%
% interpolate structure factors around reference peaks

[h,k,l] = gridarray2hkl(hklGrid);
supercell = hklGrid(1).invert.P;

LDT = proc.script.LatticeDynamicsTools('supercell',supercell,'Cell',ENM.Cell,'M',M);

sffun = LDT.calc1PSFInterpFromMap(rho,MT.Grid,MT.Basis);
[Gk,ind] = LDT.precompute1PSFs(sffun,h,k,l);

% save the structure factors
save proc/reference_calc.mat Gk ind LDT


%% functions

function hklRef = selectStrongBraggPeaks(hklTable,Basis,numRef,resMin,resMax)

% filter by reflection range
[sx,sy,sz] = Basis.invert.frac2lab(hklTable.h,hklTable.k,hklTable.l);
s = sqrt(sx.^2 + sy.^2 + sz.^2);
isIncl = 1./s >= resMin & 1./s <= resMax & ~isnan(hklTable.Fobs);
hklTable = hklTable(isIncl,:);

% sort by Fobs and choose the most intense entries
[~,ixorder] = sort(hklTable.Fobs,'descend');
hklTable = hklTable(ixorder(1:numRef),:);

hklRef = table2array(hklTable(:,{'h','k','l'}));

end

function [h,k,l] = gridarray2hkl(hklGrid)

[h,k,l] = arrayfun(@(g) g.grid(),hklGrid,'Uni',0);
h = cell2mat(cellfun(@(v) shiftdim(v,-1),h,'uni',0));
k = cell2mat(cellfun(@(v) shiftdim(v,-1),k,'uni',0));
l = cell2mat(cellfun(@(v) shiftdim(v,-1),l,'uni',0));

end
