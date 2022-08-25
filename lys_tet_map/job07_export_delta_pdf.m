%% dpdf maps for experiment

mapFileName = 'export/lys_tet_map.h5';

%% load the experimental map (variational intensity)

[MT,I] = proc.script.MapTools.import(mapFileName,'variational','I');

% symmetry expand
tmp = MT.array2table(I);
MT = MT.resize('symexpand');
I = MT.table2array(tmp,'symexpand');

clear tmp

%% calculate the dpdf

[P,MTp] = MT.fourier_transform(I);

%% export

chunksize = MTp.Grid.P; 

MTp.export('export/lys_tet_delta_pdf.h5','delta_pdf',P,'Datatype','single','ChunkSize',chunksize);

