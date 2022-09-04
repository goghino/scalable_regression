function test_bench(Xbeta_file, beta_eLR_file, Solver)

    %%% path to the IPOPT library
    addpath( ...
        'Ipopt/Ipopt/lib', ...
        '-end' );

    %%% Set the number of threads
    setenv('OMP_NUM_THREADS','1');

    %%% Some other config
    ITERATIVE = 0;
    randn('seed',1)
    MAT_WRITE = 1;

    if Solver == 3 %see config.sh for constants
        FMINCON = 1;
    else
        FMINCON = 0;
    end


    options = struct(  'dim', -1,...
                         'T', -1,...
                 'reg_param', -1,...
                  'n_repeat', -1,...
                       'tol', 1.0000e-9,...
               'mu_strategy', 'adaptive',...
                 'mu_oracle', 'probing',...
        'nlp_scaling_method', 'none',...
            'scaling_target', 1.0,...
                    'run_id', -11);
    
    Xbeta_all = dlmread(Xbeta_file);
    Dimension = Xbeta_all(1,1)
    T = Xbeta_all(2,1)
    reg_param = Xbeta_all(3,1)
    Xbeta = Xbeta_all(4:end, :);
    beta_eLR = dlmread(beta_eLR_file)';

    % NLP problem structures
    %alpha=(1/T)*X*X'.*(beta'*beta);
    alpha = Xbeta*Xbeta';

    W = randn(1,Dimension);
    fff0 = LogLik_SPACL_W(W,Dimension,alpha,beta_eLR,reg_param);

    Aeq=ones(1,Dimension);
    beq=1;

%%%%%%%%%%%%%%%%%%%%%% FMINCON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (FMINCON==1)
    fprintf("Running FMINCON\n");

    A=-speye(Dimension);b=sparse(Dimension,1);
    options = optimoptions(@fmincon,...
        'Algorithm','interior-point',...
        'MaxIter',200,...
        'SpecifyObjectiveGradient',true, ...
        'HessianFcn',@(x,lambda)hessinterior(x,lambda,alpha,reg_param),'HessPattern',0,...
        'Display','iter','TolFun',1e-20,'TolCon',1e-14,'TolPCG',1e-14,'TolX',...
        1e-14,'TolConSQP',1e-14,'TolGradCon',1e-14,'TolPCG',1e-14,'OptimalityTolerance',1e-20,'StepTolerance',1e-20);

    options.ConstraintTolerance = 1e-12;

    tic
    [W,fff11,flag11,info] =  fmincon(@(x)LogLik_SPACL_W...
        (x,Dimension,alpha,beta_eLR,reg_param)...
        ,W,(A),(b),Aeq,beq,[],[],[],options);
    toc

end
%%%%%%%%%%%%%%%%%%%%%%%%  IPOPT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (FMINCON==0)
    fprintf("Running IPOPT\n");
    W0           = W;        % The starting point.
    ipoptions.lb = zeros(1,Dimension);   % Lower bound on the variables.
    ipoptions.ub = ones(1,Dimension);    % Upper bound on the variables.
    ipoptions.cl = [0];          % Lower bounds on the constraint functions.
    ipoptions.cu = [0];          % Upper bounds on the constraint functions.


    % Set the IPOPT options.
    ipoptions.ipopt.tol             = options.tol;
    ipoptions.ipopt.dual_inf_tol    = options.tol;
    ipoptions.ipopt.constr_viol_tol = options.tol;
    ipoptions.ipopt.acceptable_tol  = options.tol;
    ipoptions.iptop.acceptable_iter = 0;
    ipoptions.ipopt.max_iter        = 200;

    ipoptions.ipopt.linear_solver      = 'pardiso';
    ipoptions.ipopt.mu_strategy        = options.mu_strategy;
    if (options.mu_strategy == "monotone")
        ipoptions.ipopt.mu_init            = options.mu_init;
    else
        ipoptions.ipopt.mu_oracle          = options.mu_oracle;
        ipoptions.ipopt.fixed_mu_oracle    = options.mu_oracle;
    end

    ipoptions.ipopt.nlp_scaling_method = options.nlp_scaling_method;
    if (options.nlp_scaling_method == "gradient-based")
        ipoptions.ipopt.nlp_scaling_max_gradient           =    options.scaling_target;
        ipoptions.ipopt.nlp_scaling_obj_target_gradient    =    options.scaling_target;
        ipoptions.ipopt.nlp_scaling_constr_target_gradient =    options.scaling_target;
    end

    % ipoptions.ipopt.corrector_type     = 'affine';
    % ipoptions.ipopt.mu_max_fact                        =    1;
    % ipoptions.ipopt.adaptive_mu_globalization          = 'obj-constr-filter';

    ipoptions.ipopt.bound_relax_factor      = 0;
    ipoptions.ipopt.honor_original_bounds   = 'yes';

    ipoptions.ipopt.jac_c_constant          = 'yes';

    ipoptions.ipopt.print_level             = 5;
    ipoptions.ipopt.print_timing_statistics = 'yes';
    ipoptions.ipopt.output_file             = 'ipopt.out';
    ipoptions.ipopt.pardiso_msglvl = 0;

    if (ITERATIVE==1)
        ipoptions.ipopt.pardiso_iterative = 'yes';
        ipoptions.ipopt.output_file       = 'ipopt-i.out';
        ipoptions.ipopt.pardiso_max_iter = 500;
        ipoptions.ipopt.pardiso_iter_relative_tol = 1e-6;
        ipoptions.ipopt.pardiso_iter_coarse_size = 5000;
        ipoptions.ipopt.pardiso_iter_max_levels = 10;
        ipoptions.ipopt.pardiso_iter_dropping_factor = 0.5;
        ipoptions.ipopt.pardiso_iter_dropping_schur = 1e-1;
        ipoptions.ipopt.pardiso_iter_max_row_fill = 10000000;
        ipoptions.ipopt.pardiso_iter_inverse_norm_factor = 5000000;
        %max number of decreases of drop tol. during one solve
        ipoptions.ipopt.pardiso_max_droptol_corrections = 4; 

        %other options:
        %DPARM_[0] = pardiso_max_iter; // Maximum number of Krylov iterations
        %DPARM_[1] = pardiso_iter_relative_tol; // Relative Residual Convergence
        %DPARM_[2] = pardiso_iter_coarse_size; // Maximum Size of Coarse Grid Matrix
        %DPARM_[3] = pardiso_iter_max_levels; // Maximum Number of Grid Levels
        %DPARM_[4] = pardiso_iter_dropping_factor;  // dropping value for incomplete factor
        %DPARM_[5] = pardiso_iter_dropping_schur;  // dropping value for sparsify schur complementfactor
        %DPARM_[6] = pardiso_iter_max_row_fill;  // max fill for each row
        %DPARM_[7] = pardiso_iter_inverse_norm_factor;  // dropping value for sparsify schur complementfactor
    end

    %IPARM(13)
    %ipoptions.ipopt.pardiso_matching_strategy = 'complete';
    %ipoptions.ipopt.pardiso_matching_strategy = 'complete+2x2';
    ipoptions.ipopt.pardiso_matching_strategy = 'constraints';

    %IPARM(2)
    %ipoptions.ipopt.pardiso_order = 'amd';
    %ipoptions.ipopt.pardiso_order = 'metis';
    ipoptions.ipopt.pardiso_order = 'pmetis';

    ipoptions.ipopt.pardiso_max_iterative_refinement_steps = 0;

    ipoptions.auxdata = struct( ...
    'd',         Dimension, ...
    'alpha',     alpha, ...
    'sp_alpha',  sparse(tril(alpha)), ...
    'beta',      beta_eLR, ...
    'eps_C',     reg_param, ...
    'Aeq',       Aeq, ...
    'beq',       beq, ...
    'MAT_WRITE', MAT_WRITE);

    % The callback functions.
    funcs.objective         = @objective;
    funcs.gradient          = @gradient;
    funcs.constraints       = @constraints;
    funcs.jacobian          = @(x, aux) sparse(Aeq);
    funcs.hessian           = @hessian;
    funcs.jacobianstructure = @(aux) sparse(Aeq);
    funcs.hessianstructure  = @(aux) sparse(tril(ones(size(W,2))));


    % Run IPOPT.
    [W info] = ipopt_auxdata(W0,funcs,ipoptions);
    %assert(info.status == 0)     
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if MAT_WRITE == 1
        fprintf('Writing data vectors into files....\n')
        FileNameNew = "W_matlab_"+num2str(Dimension)+"_"+num2str(T)+".csv";
        dlmwrite(FileNameNew,W','precision',16)
    end

end
%%
function [fun,grad]=LogLik_SPACL_W(W,d,alpha,beta,eps_C)
%fun=-2*beta*W'+W*alpha*W'+eps_C*sum(W.*(log(max(W,1e-12))));
fun=(beta+W*alpha)*W'+eps_C*sum(W.*(log(max(W,1e-12))));
if nargout > 1
    grad=real(beta+2*(alpha*W')'+eps_C.*(log(max(W,1e-12))+ones(1,d)));
end
end
 
function [grad]=LogLik_SPACL_W_grad(W,d,alpha,beta,eps_C)
grad=real(beta+2*(alpha*W')'+eps_C.*(log(max(W,1e-12))+ones(1,d)));
end

function H = hessinterior(W,lambda,alpha,eps_C)
%HMFLEQ1 Hessian-matrix product function for BROWNVV objective.
%   W = hmfleq1(Hinfo,Y,V) computes W = (Hinfo-V*V')*Y
%   where Hinfo is a sparse matrix computed by BROWNVV 
%   and V is a 2 column matrix.
H=2*alpha+diag(eps_C./max(W,1e-12));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----  callback functions  -----
function f = objective(x, aux)
%if(min(x) < 0)
%   x = max(x,0); 
%end
f = LogLik_SPACL_W(x,aux.d,aux.alpha,aux.beta,aux.eps_C);
end

function df = gradient(x, aux)
%if(min(x) < 0)
%   x = max(x,0); 
%end
df = LogLik_SPACL_W_grad(x,aux.d,aux.alpha,aux.beta,aux.eps_C);
end

function c = constraints(x, aux)
%if(min(x) < 0)
%   x = max(x,0); 
%end
c = aux.Aeq*x' - aux.beq;
end

function H = hessian(x, sigma, lambda, aux)
%if(min(x) < 0)
%   x = max(x,0); 
%end
%H1 = sparse(tril(sigma * hessinterior(x',lambda,aux.alpha,aux.eps_C)));
%H2 = tril(sigma*2*aux.alpha+diag(aux.eps_C./max(x',1e-12)));
%H = sigma*2*aux.sp_alpha+spdiags(aux.eps_C./max(x',1e-12),0,aux.d,aux.d);
H = sigma*2*aux.sp_alpha;

for i = 1:aux.d
    H(i,i) = H(i,i) + sigma*aux.eps_C./max(x(i),1e-12);
end
%H(1:(aux.d+1):end) = full(H(1:(aux.d+1):end))' + aux.eps_C./max(x',1e-12);

end
