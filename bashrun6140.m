function [] = bashrun6140(i)

pat = '6140'; %patient number

raw = load(strcat('pat',pat,'.mat'));
cell = struct2cell(raw);
K = size(cell,1); %no. of variables
N = size(cell{1},2); %no. of observations
M = zeros(N,K);
for j = 1:K
    M(:,j) = cell{j};
end

%extract the bit to use
Mat = M(:,:);
n = size(Mat,1);
dims = size(Mat,2);

[RM,ecC] = PMIMEsig(Mat(i+1:i+1000,:),5,1,5,100,0.05,0);
dlmwrite(strcat('P_mat_pat',pat,'_',int2str(i),'.txt'),RM,'delimiter','\t');

end
