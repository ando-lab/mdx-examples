%% generate report

reportDir = 'report';
reportFileName = 'report_lys_tri_model.md';
mkdir(reportDir)

%%
% initialize markdown file

report_md = {'# Report: triclinic lysozyme modeling',datestr(now)};

%% 

load proc/reference_halos.mat I sigma
fh = figure(1);clf;fh.Position(3:4) = [425,320];

ndata = nnz(~isinf(sigma));
I(isinf(sigma)) = NaN;

m = cell(size(I,1),1);
for j=1:size(I,1)
    m{j} = squeeze(I(j,:,:,ceil(size(I,4)/2)));
end
imagesc(cell2mat(reshape(m,16,[])),[0,10*std(I(~isnan(I)))]);
set(gca,'Position',[0,0,1,1]);
axis image;axis off;colormap turbo

figName = 'reference_halos_lys_tri.png';
exportgraphics(fh,fullfile(reportDir,figName));

report_md = [report_md,{'## Reference halos'}];
report_md = [report_md,{sprintf('![](%s)',figName)}];
%%
% plot GOODVIBES refinement info

load proc/goodvibes_fit.mat refineLog
load proc/reference_halos.mat sigma
ndata = nnz(~isinf(sigma));

md_table = {'Stage | Spring type | Parameterization | Number of parameters | Number of iterations | Chi-squared',...
'--- | --- | --- | --- | --- | ---'};

for j=1:numel(refineLog)
    st = refineLog(j).springType;
    pa = refineLog(j).parameterization;
    np = numel(refineLog(j).pfit);
    ni = refineLog(j).history(end).iteration;
    x2 = refineLog(j).history(end).resnorm/ndata;
    md_table = [md_table,{sprintf('%d | %s | %s | %d | %d | %g',j,st,pa,np,ni,x2)}];
end

report_md = [report_md,{'## GOODVIBES refinement',strjoin(md_table,'\n')}];

%% compute intensity statistics

slim = 0.8;
edges = 0.04:.01:slim;

mapFileName = '../lys_tri_map/export/lys_tri_map.h5';

[MT,I,I1,I2] = proc.script.MapTools.import(mapFileName,'total','I','I1','I2');

Ilatt = MT.read_data('export/lys_tri_goodvibes.h5','supercell','I');

T = calc_intensity_stats(MT,I,I1,I2,Ilatt,edges);

%% plot intensity statistics

fh = figure(1);clf;fh.Position(3:4) = [300,300];
[t,ah] = plot_intensity_stats(T);

figName = 'intensity_stats_lys_tri.png';
exportgraphics(fh,fullfile(reportDir,figName));

report_md = [report_md,{'## Intensity statistics'}];
report_md = [report_md,{sprintf('![](%s)',figName)}];

%% GOODVIBES joint ADPs

Tdb = readtable('export/lys_tri_goodvibes_com_covariances.csv');
T = preprocess_covariance_table(Tdb);

fh = figure(1);clf;fh.Position(3:4) = [350,250];
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;scatter(T.d,T.viso,5,'Filled','MarkerFaceColor','b');
ylabel('Total (Å^2)');box on;set(gca,'XTickLabel',{});

nexttile;scatter(T.d,T.vaniso,5,'Filled');
lg = legend({'1,1','2,2','3,3','1,2','1,3','2,3'});
lg.Location = 'bestoutside';
xlabel('distance (Å)');ylabel('Anisotropic (Å^2)');box on;

figName = 'goodvibes_joint_adps_lys_tri.png';
exportgraphics(fh,fullfile(reportDir,figName));

report_md = [report_md,{'## GOODVIBES joint-ADPs'}];
report_md = [report_md,{sprintf('![](%s)',figName)}];
%% DISCOBALL analysis (iso vs. distance)

% iso covariance vs distance
Tdb = readtable('export/lys_tri_discoball_joint_adps.csv');
T = preprocess_covariance_table(Tdb);

fh = figure(1);clf;fh.Position(3:4) = [350,250];
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;scatter(T.d,T.viso,5,'Filled','MarkerFaceColor','b');
ylabel('Total (Å^2)');box on;set(gca,'XTickLabel',{});

nexttile;scatter(T.d,T.vaniso,5,'Filled');
lg = legend({'1,1','2,2','3,3','1,2','1,3','2,3'});
lg.Location = 'bestoutside';
xlabel('distance (Å)');ylabel('Anisotropic (Å^2)');box on;

figName = 'discoball_joint_adps_lys_tri.png';
exportgraphics(fh,fullfile(reportDir,figName));

report_md = [report_md,{'## DISCOBALL joint-ADPs'}];
report_md = [report_md,{sprintf('![](%s)',figName)}];


%% DISCOBALL validation

Tdb = readtable('export/lys_tri_discoball_joint_adps.csv');
Tgv = readtable('export/lys_tri_goodvibes_joint_adps.csv');

[T,fit_iso,fit_aniso] = compare_covariance_tables(Tdb,Tgv);

fh = figure(1);clf;fh.Position(3:4) = [350,180];
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;scatter_goodvibes_vs_discoball(fit_iso);title('isotropic');
nexttile;scatter_goodvibes_vs_discoball(fit_aniso);title('anisotropic');

figName = 'discoball_validation_lys_tri.png';
exportgraphics(fh,fullfile(reportDir,figName));

report_md = [report_md,{'## DISCOBALL validation'}];
report_md = [report_md,{sprintf('![](%s)',figName)}];

%% B-factors

load('proc/atomic_model.mat','AtomFF','Atoms');
Atoms.mdxFormFactor = AtomFF.mdxFormFactor;
load('proc/goodvibes_fit.mat','LD');
load('proc/goodvibes_model.mat','ENM');
tlsori = ENM.Cell.AsymmetricUnit.ori;

C = LD.cov_unitcell;
[T,L,S,rcor,Tcor,Lcor,Scor] = tls_analysis(C);

md_matrix = @(m) strjoin([{'```'},cellstr(num2str(m))',{'```'}],'\n');

report_md = [report_md,{'## GOODVIBES ADPs'}];
report_md = [report_md,{...
    '**Center of mass** (Å)',md_matrix(tlsori(:)'),...
    '**Center of reaction** (Å)',md_matrix(rcor(:)' + tlsori(:)'),...
    '**T** (center of reaction, Å^2)',md_matrix(Tcor),...
    '**L** (center of reaction, deg^2)',md_matrix(Lcor*(180/pi)^2),...
    '**S** (center of reaction, Å deg)',md_matrix(Scor*(180/pi))}];

stats = calc_backbone_bfactors(Atoms,C,tlsori,Tcor);

fh = figure(1);clf;fh.Position(3:4) = [420,220];
t = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
nexttile;
[ah,a] = Bplot_stacked(stats);
xlabel('Residue number');
ylabel('Mean B-factor of backbone atoms (Å^2)');
lg = legend(a,{'Total','Residual','GOODVIBES: Rotation','GOODVIBES: Translation'},'Location','northwest','Box','off');

figName = 'bfactors_lys_tri.png';
exportgraphics(fh,fullfile(reportDir,figName));

report_md = [report_md,{sprintf('![](%s)',figName)}];

%% delta pdf slices

rmin = 2.5;
rmax = 25; % max radius

mapFileName = 'export/lys_tri_goodvibes_delta_pdf.h5';

MT = proc.script.MapTools.import(mapFileName,'delta_pdf');
MT = MT.resize('radius',rmax);
[MTx,MTy,MTz] = orthoslice_grid(MT);

Pz_sim = MTz.read_data(mapFileName,'delta_pdf','P');

mapFileName = '../lys_tri_map/export/lys_tri_delta_pdf.h5';
Pz_exp = MTz.read_data(mapFileName,'delta_pdf','P');

z = numel(MT.SpaceGroup.generalPositions);
Pz_sim = Pz_sim/z;
Pz_exp = Pz_exp/z;

msk = MTz.spherical_mask(rmin) | ~MTz.spherical_mask(rmax);
Pz_exp(msk) = 0;
Pz_sim(msk) = 0;

clim = [-50,50];
fh = figure(1);clf;fh.Position(3:4) = [400,145];
t = tiledlayout(1,3,'TileSpacing','tight','Padding','tight');
nexttile;imagesc(squeeze(Pz_exp)',clim);colormap gray;axis image;axis off;title('Experiment')
nexttile;imagesc(squeeze(Pz_sim)',clim);colormap gray;axis image;axis off;title('GOODVIBES')
nexttile;imagesc(squeeze(Pz_exp-Pz_sim)',clim);colormap gray;axis image;axis off;title('residual');colorbar
get(gca,'Clim')


figName = 'delta_pdf_subtraction.png';
exportgraphics(fh,fullfile(reportDir,figName));

report_md = [report_md,{'## 3D-ΔPDF subtraction'}];
report_md = [report_md,{sprintf('![](%s)',figName)}];

%% delta pdf stats

rmin = 2.5;
rmax = 25;
dr = 1;

mapFileName = 'export/lys_tri_goodvibes_delta_pdf.h5';

MT = proc.script.MapTools.import(mapFileName,'delta_pdf');
MT = MT.resize('radius',rmax);

nsymop = numel(MT.SpaceGroup.generalPositions);

P_sim = MT.read_data(mapFileName,'delta_pdf','P');
P_sim = P_sim/nsymop; % per asymmetric unit

mapFileName = '../lys_tri_map/export/lys_tri_delta_pdf.h5';
P_exp = MT.read_data(mapFileName,'delta_pdf','P');
P_exp = P_exp/nsymop; % per asymmetric unit

stats = calc_dpdf_stats(MT,P_exp,P_sim,rmin:dr:rmax);

fh = figure(1);clf; fh.Position(3:4) = [250,180];
plot_delta_pdf_stdev(stats);

figName = 'delta_pdf_statistics.png';
exportgraphics(fh,fullfile(reportDir,figName));

report_md = [report_md,{'## 3D-ΔPDF statistics'}];
report_md = [report_md,{sprintf('![](%s)',figName)}];

%% Diffuse scattering subtraction


mapFileName = 'export/lys_tri_map_sub.h5';

MT = proc.script.MapTools.import(mapFileName,'sub_fill');
[MTx,MTy,MTz] = orthoslice_grid(MT);

[I_sub_fill,sigma_fill] = MTz.read_data(mapFileName,'sub_fill','I','sigma');
I_sub_fill = MTz.table2array(MTz.array2table(I_sub_fill),'symexpand','replace',NaN);

mapFileName = '../lys_tri_map/export/lys_tri_map.h5';
[I_exp] = MTz.read_data(mapFileName,'total','I');
I_exp = MTz.table2array(MTz.array2table(I_exp),'symexpand','replace',NaN);

mapFileName = 'export/lys_tri_goodvibes.h5';
[I_sim] = MTz.read_data(mapFileName,'supercell','I');
I_sim = MTz.table2array(MTz.array2table(I_sim),'symexpand','replace',NaN);

z = numel(MT.SpaceGroup.generalPositions);
I_sub_fill = I_sub_fill/z;
I_exp = I_exp/z;
I_sim = I_sim/z;

clim = [0,0.8E5];
fh = figure(1);clf;fh.Position(3:4) = [400,400];
t = tiledlayout(2,2,'TileSpacing','tight','Padding','tight');
nexttile;imagesc(squeeze(I_exp)',clim);colormap turbo;axis image;axis off;title('Experiment')
nexttile;imagesc(squeeze(I_sim)',clim);colormap turbo;axis image;axis off;title('GOODVIBES')
nexttile;imagesc(squeeze(I_exp - I_sim)',clim);colormap turbo;axis image;axis off;title('residual');
nexttile;imagesc(squeeze(I_sub_fill)',clim);colormap turbo;axis image;axis off;title('filtered residual');colorbar
get(gca,'Clim')

figName = 'diffuse_subtraction.png';
exportgraphics(fh,fullfile(reportDir,figName),'Resolution',300);

report_md = [report_md,{'## Diffuse subtraction'}];
report_md = [report_md,{sprintf('![](%s)',figName)}];

%% write report file
fid = fopen(fullfile(reportDir,reportFileName),'w');
fprintf(fid,'%s\n\n',report_md{:});
fclose(fid);
%% FUNCTIONS

function [MTx,MTy,MTz] = orthoslice_grid(MT)

% prepare slices
[f1,f2,f3] = MT.Grid.ind2frac(MT.Grid.N(1),MT.Grid.N(2),MT.Grid.N(3));
[MTx] = MT.resize('roi',[0,0,-f2,f2,-f3,f3]);
[MTy] = MT.resize('roi',[-f1,f1,0,0,-f3,f3]);
[MTz] = MT.resize('roi',[-f1,f1,-f2,f2,0,0]);

end

function plot_delta_pdf_stdev(stats)
c = lines(7);
rmax = max(stats.distance) + 0.5*diff(stats.distance(1:2));
scatter(stats.distance,stats.total,25,'ko');hold on;
scatter(stats.distance,stats.goodvibes,25,'<','Filled','MarkerFaceColor',c(1,:));
scatter(stats.distance,stats.residual,25,'d','Filled','MarkerFaceColor',c(5,:));

set(gca,'Xlim',[0,rmax],'Xscale','Lin','Yscale','Log');

xlabel('Pair distance (Å)')
ylabel({'Standard deviation of','3D-\DeltaPDF (per ASU)'})
legend({'total','GOODVIBES','residual'},'Box','off')

box on;
end



function dpdf_stats = calc_dpdf_stats(MT,dpdf_exp,dpdf_sim,edges)

% compute radial stats
SVR = proc.script.StatisticsVsRadius(...
    'edges',edges,...
    'Basis',MT.Basis,...
    'PeriodicGrid',MT.Grid);

Tsim = SVR.run(dpdf_sim.^2);
Texp = SVR.run(dpdf_exp.^2);
Tresid = SVR.run((dpdf_exp-dpdf_sim).^2);

dpdf_stats = table(Tsim.r,sqrt(Tsim.av),sqrt(Texp.av),sqrt(Tresid.av),...
    'VariableNames',{'distance','goodvibes','total','residual'});

end

function scatter_goodvibes_vs_discoball(stats)

xmin = min([min(stats.x),min(stats.y)]);
xmax = max([max(stats.x),max(stats.y)]);

xlims = [xmin,xmax]*1.05;
scatter(stats.x,stats.y,5,'Filled','MarkerFaceColor','b');
hold on;plot(xlims,xlims,'-','Color',[.6,.6,.6]);
hold on;plot(xlims,polyval(stats.p,xlims),'-','Color','k')
axis equal;
box on;
set(gca,'Xlim',xlims,'Ylim',xlims,'TickLength',[0.02,0]);
text(gca,.1,.9,sprintf('r = %.4f',stats.cc),'Units','normalized','Color','k');
text(gca,.1,.8,'y = x','Units','normalized','Color',[.6,.6,.6]);

ylabel('GOODVIBES');
xlabel('DISCOBALL');

end

function T = calc_intensity_stats(MT,I,I1,I2,I_latt,edges)

SVR = proc.script.StatisticsVsRadius(...
    'edges',edges,...
    'Basis',MT.Basis,...
    'PeriodicGrid',MT.Grid);

stats_exp_half = SVR.run(I1,I2);
stats_ld_exp = SVR.run(I_latt,I);

cc12 = stats_exp_half.cc;
cc12(cc12<0) = 0; % make sure ccstar is real
ccstar = sqrt(2*cc12./(1+cc12));

s = stats_ld_exp.r;
sigmaI_exp = sqrt(stats_ld_exp.var2).*ccstar;
sigmaI_latt = sqrt(stats_ld_exp.var1);
cclatt = stats_ld_exp.cc;

T = table(s,sigmaI_exp,sigmaI_latt,cclatt,ccstar);

end

function [t,ah] = plot_intensity_stats(T)

t = tiledlayout(2,1);
%t.XLabel.String = '1/d (Å^{-1})';
t.TileSpacing='tight';
t.Padding='compact';

bin_width = T.s(end)-T.s(end-1);
slim = T.s(end) + bin_width/2;
ymax = 1.2*max(T.sigmaI_exp);

nexttile;
plot(T.s,T.sigmaI_exp,'k.-')
hold on;
plot(T.s,T.sigmaI_latt,'m.-')
ah(1) = gca;
set(ah(1),'Ylim',[0,ymax],'Xlim',[0,slim],'XTickLabel',{},'Box','Off','TickDir','out')
legend({'expt.','sim.'},'Location','northeast','Box','Off')
ylabel('Standard deviation')

nexttile;
plot(T.s,T.ccstar,'k-','Color',[1,1,1]*.5)
hold on;
plot(T.s,T.cclatt,'k.-');
plot([0,slim],[1,1],'k--');
ah(2) = gca;
set(ah(2),'Ylim',[0.5,1.05],'Xlim',[0,slim],'Box','Off','TickDir','out')
legend({'CC*','expt. vs. sim.'},'Location','south','Box','Off')
ylabel('Correlation coefficient')
xlabel('1/d (Å^{-1})')

end

function T = preprocess_covariance_table(T)
T.d = sqrt(T.x.^2 + T.y.^2 + T.z.^2);
isIncl = T.d>0;
T = T(isIncl,:);
T.viso = T.v11 + T.v22 + T.v33;
T.vaniso = cat(2,T.v11 - T.viso/3,T.v22 - T.viso/3, T.v33 - T.viso/3, T.v12,T.v13,T.v23);
T = T(:,{'n1','n2','n3','d','viso','vaniso'});
end

function [T,fit_iso,fit_aniso] = compare_covariance_tables(Tdb,Tgv)

discoball = preprocess_covariance_table(Tdb);
goodvibes = preprocess_covariance_table(Tgv);
T = innerjoin(discoball,goodvibes,'Keys',{'n1','n2','n3'},'RightVariables',{'viso','vaniso'});

x = T.viso_discoball;
y = T.viso_goodvibes;

p = polyfit(x,y,1);
ccmat = corrcoef(x,y);
cc = ccmat(1,2);
fit_iso = struct('x',x,'y',y,'cc',cc,'p',p);

x = T.vaniso_discoball(:);
y = T.vaniso_goodvibes(:);

p = polyfit(x,y,1);
ccmat = corrcoef(x,y);
cc = ccmat(1,2);
fit_aniso = struct('x',x,'y',y,'cc',cc,'p',p);

end

function [T,L,S,rcor,Tcor,Lcor,Scor] = tls_analysis(C)

C = C(1:6,1:6);

T = C(1:3,1:3); % A^2
L = C(4:6,4:6);
S = C(4:6,1:3);

A = trace(L)*eye(3)-L;
b = [S(2,3)-S(3,2);S(3,1)-S(1,3);S(1,2)-S(2,1)];
rcor = A\b;

cmat = @(r1,r2,r3) sparse([3,1,2,2,3,1],[2,3,1,3,1,2],[r1,r2,r3,-r1,-r2,-r3],3,3);

Oshift = kron(sparse(1,2,1,2,2),cmat(-rcor(1),-rcor(2),-rcor(3))) + speye(6,6);

C0 = Oshift*C*Oshift';

Tcor = C0(1:3,1:3); % A^2
Lcor = C0(4:6,4:6);
Scor = C0(4:6,1:3);

end

function [stats] = calc_backbone_bfactors(Atoms,C,tlsori,Tcor)

% per-residue B-factors (BACKBONE atoms)
isIncl = Atoms.mdxChemicalGroup=="backbone" & Atoms.mdxAtomicSymbol~="H";
Atoms = Atoms(isIncl,:);

[Atoms.Biso] = calc_Biso(Atoms.mdxFormFactor);
[Atoms.Blatt] = calc_Blatt(Atoms,C,tlsori); 

stats = groupsummary(Atoms,{'chainID','resSeq'},{'min','max','mean'},...
    {'Biso','Blatt'},'IncludeEmptyGroups',true);

stats.Blatt_cor = ones(size(stats,1),1)*8*pi^2*trace(Tcor)/3;

function B = calc_Biso(FF)
U = arrayfun(@(ff) ff.U,FF,'Uni',0);
B = cellfun(@(u) 8*pi^2*trace(u)/3,U);
end

function Blatt = calc_Blatt(Atoms,Clatt,ori)
if nargin < 3 || isempty(ori)
    ori = [0,0,0];
end
P = arrayfun(@(x,y,z) nm.Group(x,y,z).tl2uxyz,Atoms.x - ori(1),Atoms.y - ori(2),Atoms.z - ori(3),'Uni',0);
Ulatt = cellfun(@(p) p*Clatt(1:6,1:6)*p',P,'Uni',0);
Blatt = cellfun(@(u) 8*pi^2*trace(u)/3,Ulatt);
end


end

function [ah,a] = Bplot_stacked(stats)

% colors: https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=3
c1 = [166,206,227]/255;
c2 = [31,120,180]/255;
c3 = [178,223,138]/255;

% data to plot
x = stats.resSeq;
ylatt = stats.mean_Blatt;
ylattcor = stats.Blatt_cor;
yobs = stats.mean_Biso;

% make the plot
a(1) = area(x,yobs); hold on;
a(2) = area(x,ylatt);
a(3) = area(x,ylattcor);

p = plot(x,yobs,'k.-');

% set the colors
a(1).FaceColor = c3;
a(1).EdgeColor = [0,0,0];
a(1).LineWidth = 0.5;
a(2).FaceColor = c2;
a(2).EdgeColor = [0,0,0];
a(3).FaceColor = c1;
a(3).EdgeColor = [0,0,0];

ah = gca;

% set the axis limit
set(ah,'Ylim',[0,Inf],'Xlim',[0,Inf]);
a = [p,a];

end
