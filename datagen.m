pats = {'6140', '6227', '6232', '6255', '6383a', '6383b', '6395', '6396a', '6396b', '6527', '7063', '7574', '7577', '7608', '7634', '7771', '7890', '7943'};
ns = {284, 300, 275, 330, 173, 173, 73, 268, 284, 258, 448, 295, 320, 364, 176, 290, 293, 291};
unthresholded = 'TRUE';
if unthresholded == 'TRUE'
    cd Unthresholded
end

for i = 1:length(ns) %for each patient
    pat = pats{i};   %set patient number
    n = ns{i}*1000;  %set max matrix number
    for batch = 0%:1
        dat = [];         %initialise matrix to contain metric values
        if batch==0
            letter = 'T'; %for tigramite matrices
        else
            letter = 'P'; %for PMIME matrices
        end
        for val = 0:1000:n %for each matrix
            fprintf('Progress: %0.1f%% (patient %s, batch %s) \n', val/n*100, pat, letter)
            try
                A = importdata(strcat(letter,'mat_pat',pat,'_',int2str(val),'.txt')); %import matrix
            catch %if import fails:
                fprintf('Matrix %d not found for patient %s (batch %s) \n', val, pat, letter);
                continue; %skip this loop iteration
            end

            
            L = zeros(size(A));                          %initialise connection LENGTH matrix
            ind=find(A); L(ind)=1./A(ind);               %specify lengths by inverting weights 


            [in_deg,out_deg,deg] = degrees_dir(A);       %node degrees
            in_deg_mean = mean(in_deg);
            in_deg_std = std(in_deg);
            out_deg_mean = mean(out_deg);
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
            out_str_mean = mean(out_str);
            out_str_std = std(out_str);

            str_dif = in_str - out_str;                  %difference between in and out strength
            str_dif_std = std(str_dif);

            [EBC,NBC] = edge_betweenness_wei(L);         %edge and node betweenness
            EBC_mean = mean(EBC, 'all');
            EBC_std = std2(EBC);
            NBC_mean = mean(NBC);
            NBC_std = std(NBC);

            dens = density_dir(A);                       %density i.e. wiring cost if ignoring weights
            K = sum(sum(L));                             %wiring cost = sum of connection lengths
            K = K/nnz(L);                                %normalise by number of connections

            E_glob = efficiency_wei(A, 0);               %global efficiency
            E_loc = efficiency_wei(A, 2);                %local (node) efficiencies
            E_loc_mean = mean(E_loc);
            E_loc_std = std(E_loc);

            [comms,Q] = modularity_dir(A);               %optimal community structure and modularity

            C = clustering_coef_wd(A);                   %node clustering coefficient
            C_mean = mean(C);
            C_std = std(C);

            R_oi = assortativity_wei(A,1);               %correlation of out vs in strength
            R_io = assortativity_wei(A,2);               %correlation of in vs out strength
            R_oo = assortativity_wei(A,3);               %correlation of out vs out strength
            R_ii = assortativity_wei(A,4);               %correlation of in vs in strength

            %E_cost = E_glob - K;                        %cost efficiency (removed as we create later) 

            vars = [in_deg_mean, in_deg_std, out_deg_mean, out_deg_std, deg_dif_std, in_str_mean, in_str_std, out_str_mean, out_str_std, str_dif_std, EBC_mean, EBC_std, NBC_mean, NBC_std, dens, K, E_glob, E_loc_mean, E_loc_std, Q, C_mean, C_std, R_oi, R_io, R_oo, R_ii];
            dat = [dat;vars];
        end
        dlmwrite(strcat(letter,'metrics_pat',pat,'.txt'),dat,'delimiter','\t');
    end
end
