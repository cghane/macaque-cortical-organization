%% computes voxel-wise functional connectivity matrices
%monkeys = {'Andrea', 'Daphne', 'Bilbo', 'Aragorn', 'Tarwi'};
base_path = '/mnt/scratch/NHP4CYRUS/data_dump/';
base_pathout = '/mnt/scratch/NHP4CYRUS/';
folder_info = dir(base_path);
folder_names = {folder_info([folder_info.isdir]).name}; % Only directories
monkeys = folder_names(~ismember(folder_names, {'.', '..'})); % Exclude '.' and '..'
%base_path = '/users/cyrusghane/MATLAB/Projects/FC-INT/';
% PATH_I = '/mnt/scratch/NHP4ANA/ACC/fMRI/mat_files/';
%PATH_I = '/users/cyrusghane/MATLAB/Projects/FC-INT/Tarwi';
%PATH_O = '/users/cyrusghane/MATLAB/Projects/FC-INT/FCdata/FCvisual/';
%PATH_O = '/mnt/data/NHP4ANA/striatum_frontal/';
% Import cort_crd and subcort_crd
% These are the coordinates of the hippocampus for example
area_name = 'cortical';
for i = 1:length(monkeys)
    monkey = monkeys{i};
    PATH_I = [base_path, monkey];
    PATH_O = [base_pathout, 'FCdata/FCVisual/'];
    conns = {};
    % Import data
    for run = 1:5
        data_cortical = importdata(['/mnt/scratch/NHP4CYRUS/data_dump/', monkey, '/', ...
                                   area_name, '/ROIdump_', area_name, num2str(run), '.1D']);
                 %disp(['Loading: ', char(PATH_I), '/subcortical/ROIdump_',...
                                      %'subcortical' , num2str(run), '.1D']); 
        %data_subcortical = importdata([char(PATH_I), '/subcortical/ROIdump_',...
                                  % 'subcortical' , num2str(run), '.1D']);

        %datacat_cortical{run} = data_cortical(crd_acc,:);
        datacat_cortical{run} = data_cortical;
        %datacat_subcortical{run} = data_subcortical(crd_acc,:);
    end
    clearvars data_cortical data_subcortical
    for k = 1:5
      %conn = corr(datacat_subcortical{k}', datacat_cortical{k}');     
      %conn = corr(datacat_subcortical{k}');
      conn = corr(datacat_cortical{k}');
      conns{k} = single(atanh(conn));
    end
    
    %save([PATH_O, 'andreahippocampus','.mat'], 'conns', 'cort_crd','subcort_crd','-v7.3')
    save([PATH_O, lower(monkey),'.mat'], 'conns','crd_acc','-v7.3')
    conns = []; datacat_cortical =[]; datacat_subcortical =[];
    %average connectivity matrices across runs%
end
