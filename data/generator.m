MAT_WRITE = 1;
T=1000;
Dimension=[1024 2048 4096 8192 16384 2*16384 4*16384 8*16384 16*16384];
reg_param=1e-3;

randn('seed',1)
for d=1:length(Dimension)
    disp(['dimension=' num2str(Dimension(d))]);
    
    % Inputs
    X=randn(Dimension(d),T);
    pi=randn(1,T);
    beta=randn(1,Dimension(d));
     
    % Auxiliary structure
    Xbeta = X;
    for i=1:T
        Xbeta(:,i) = X(:,i).*beta'*sqrt(1/T);
    end

    % NLP problem structures
    %alpha=(1/T)*X*X'.*(beta'*beta);
    alpha = Xbeta*Xbeta';
    e = ones(T,1);
    beta_eLR=(-2/T)*beta.*((pi.*X) * e)';

    if MAT_WRITE == 1
        %%Xbeta
        FileName = 'Xbeta.csv';
        dlmwrite(FileName,Xbeta,'precision',16);
        
        %%stats, gradient, Xbeta
        S = fileread(FileName);
        S = [num2str(Dimension(d)), newline, num2str(T), newline, num2str(reg_param), newline, S];
        FileNameNew = "Xbeta_"+num2str(Dimension(d))+"_"+num2str(T)+".csv";
        FID = fopen(FileNameNew, 'w');
        if FID == -1, error('Cannot open file %s', FileNameNew); end
        fwrite(FID, S, 'char');
        fclose(FID);
        
        %%beta_eLR
        FileNameNew = "beta_eLR_"+num2str(Dimension(d))+"_"+num2str(T)+".csv";
        dlmwrite(FileNameNew,beta_eLR','precision',16);
    end
end
