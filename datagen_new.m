%pats = {'6140', '6227', '6232', '6255', '6383a', '6383b', '6395', '6396a', '6396b', '6527', '7063', '7574', '7577', '7608', '7634', '7771', '7890', '7943'};
%ns = {284, 300, 275, 330, 173, 173, 73, 268, 284, 258, 448, 295, 320, 364, 176, 290, 293, 291}; %num of readings per patient/1000
%unthresholded = 'TRUE';
%if unthresholded == 'TRUE'
%    cd Unthresholded
%end
BCTFolder=strcat(pwd,'\2019_03_03_BCT');
addpath(BCTFolder);
schiz = [1461, 2117, 3464, 3757, 6527, 6568, 7063, 7574, 7608, 7634, 7771, 7943];
health = [2261, 2645, 6140, 6227, 6232, 6255, 6383, 6395, 6396, 7577, 7890, 8368];
source = strcat(pwd,'\Networks\');
data_files = dir(strcat(source,'*.txt'));
F = length(data_files);
%fprintf(int2str(F));
dat = zeros([F,25]);
for i = 1:F %for each patient
    fprintf(strcat(int2str(i),' out of\t',int2str(F),'\n'));
    fname = data_files(i).name;
    parts = split(fname,'_');
    id = str2num(parts{1});
    if strcmp(parts{2},'RestEyesClosed')
        state = 0;
    else
        state = 1;
    end
    if any(schiz(:)==id)
        sch = 1;
    else
        sch = 0;
    end
    A = importdata(strcat(source,fname));
    A = A(1:47,1:47);
    L = zeros(size(A));                          %initialise connection LENGTH matrix
    ind=find(A); L(ind)=1./A(ind);               %specify lengths by inverting weights
    
    
    [in_deg,out_deg,deg] = degrees_dir(A);       %node degrees
    in_deg_mean = mean(in_deg);
    in_deg_std = std(in_deg);
    %out_deg_mean = mean(out_deg);
    out_deg_std = std(out_deg);
    
    deg_dif = in_deg - out_deg;                  %difference between in and out degree
    deg_dif_std = std(deg_dif);
    
    [in_str,out_str,str] = strengths_dir(A);     %node strengths
    in_str = in_str./in_deg;                     %normalise by degree
    in_str(isnan(in_str))=0;                     %set 0/0 to 0 instead of NaN
    out_str = out_str./out_deg;
    out_str(isnan(out_str))=0;
    
    in_str_mean = mean(in_str);
    in_str_std = std(in_str);
    %out_str_mean = mean(out_str);
    out_str_std = std(out_str);
    
    str_dif = in_str - out_str;                  %difference between in and out strength
    str_dif_std = std(str_dif);
    
    [EBC,NBC] = edge_betweenness_wei(L);         %edge and node betweenness
    EBC_mean = mean2(EBC);%, 'all');
    EBC_std = std2(EBC);
    NBC_mean = mean(NBC);
    NBC_std = std(NBC);
    
    dens = density_dir(A);                       %density i.e. wiring cost if ignoring weights
    K = sum(sum(L));                             %wiring cost = sum of connection lengths
    K = K/nnz(L);                                %normalise by number of connections
    
    E_glob = efficiency_wei(A, 0);               %global efficiency
    E_loc = efficiency_wei(A, 2);                %local (node) efficiencies
    %E_loc_mean = mean(E_loc);
    %E_loc_std = std(E_loc);
    
    [comms,Q] = modularity_dir(A);               %optimal community structure and modularity
    
    C = clustering_coef_wd(A);                   %node clustering coefficient
    C_mean = mean(C);
    C_std = std(C);
    
    R_oi = assortativity_wei(A,1);               %correlation of out vs in strength
    R_io = assortativity_wei(A,2);               %correlation of in vs out strength
    R_oo = assortativity_wei(A,3);               %correlation of out vs out strength
    R_ii = assortativity_wei(A,4);               %correlation of in vs in strength
    
    %E_cost = E_glob - K;                        %cost efficiency (removed as we create later)
    vars = [id,state,in_deg_mean,in_deg_std,out_deg_std,deg_dif_std,in_str_mean,in_str_std,out_str_std,str_dif_std,EBC_mean,EBC_std,NBC_mean,NBC_std,dens,K,E_glob,Q,C_mean,C_std,R_oi,R_io,R_oo,R_ii,sch];
    vars = real(vars);
    dat(i,:) = vars;
end
dlmwrite('metricsdata_pc_recon.txt',dat,'delimiter','\t');
