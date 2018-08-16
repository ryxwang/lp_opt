function [data_sparse_list, start_ind] = read_data(hic_files, bin_size)

%read in hi-c matrices
L=length(hic_files);
for i=1:L
    data = dlmread(hic_files{i});
    data(:,1) = data(:,1)/bin_size+1;
    data(:,2) = data(:,2)/bin_size+1;
    n_min{i}=min(data(:,1));
    n_max{i}=max(data(:,2));
    data_list{i}=data;
end
n_min=min(n_min{:}); start_ind=n_min;
n_max=max(n_max{:});
n=n_max-n_min+1;

for i=1:L
    data_list{i}(:,1) = data_list{i}(:,1)-n_min+1;
    data_list{i}(:,2) = data_list{i}(:,2)-n_min+1;
    data_sparse = sparse(data_list{i}(:,1), data_list{i}(:,2), data_list{i}(:,3), n, n);
    data_sparse = data_sparse+data_sparse'-diag(diag(data_sparse));
    data_sparse_list{i} = data_sparse;
end


%symmetry and zero diagonal




