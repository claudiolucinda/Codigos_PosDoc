clear;

% 1. Simulation settings
% Vector of seeds that creates datasets
R           =  1000;     % Number of MC simulations
seedvec     = (1:100000/R:100000);

% 1.2. TRUE PARAMETER VALUES
% Size of the datasets
nmkt        = 50;                               % number of markets = (# of cities)*(# of quarters)
prods       = 25;                               % number of brands per market
nseg        = 2;                                % number of segments in each market
nobs        = nmkt*prods;                       % number of observations
meanc       = [0 0];                            % mean of the generated vectors:
% the first is the segment latent variable, the second is the product attribute

%% Settings

% True model parameter values
betatrue    = [-1 ; -2 ; -3];                   % true mean tastes on constant, X1, dummy for segment
rc_true     = 1 ;                             % true random coefficient
SigsegTrue  = 0.3 ;

%% Set correlation: parameter that determines whether the segment is informative about the continuous characteristic or not
corrDstarX1  = 0;

%% This changes the numerosity of the segments
cutoff       = 0;

stdevDstar  = 1;                                % std. dev. of dstar (does not matter)
stdevX1     = 1;                                % std. dev. of X1

Chelp1      = [stdevDstar^2, corrDstarX1*stdevX1*stdevDstar; corrDstarX1*stdevX1*stdevDstar, stdevX1^2];    %varcovar matrix

rng default
ksi     = randn(nobs,1);
ksi     = ksi-mean(ksi);                        % rescale to have E[ksi]=0
const   = ones(nobs,1);

NonLinTheta = [SigsegTrue;rc_true];

% Number of parameters
nlin        = size(betatrue,1);                 % number of linear parameters
nrc         = size(rc_true,1);                  % number of random coeficients - Random coefficient Logit
nrcNL       = size(NonLinTheta,1);              % number of random coeficients - Random Coefficient Nested Logit

% Individual taste,integration of market shares using Quadrature Rule
% Note: when using quadrature we have no variation across markets
% In other words there is no randomness when integrating marketshares we
% use fixed nodes and weights as given by Heiss and Winschel (2007)
[nodedraws qweight] = nwspgr('KPN', nrc, 7);
Nnodes              = length(qweight);
qweight             = qweight';
qweightrprods       = reshape((qweight(ones(prods,1),:)),prods,1,1,9);

% 2. Start Simulation

% Store results
% The order of the coefficients is:
%[Constant;Segment Coefficient; Mean value of X1 (Product Attribute); Nesting parameter; Random Coefficient]
% 5=number of coefficients; 1=number of columns; 4=number of estimators (Log,NL,RC,RCNL);R=number of draws

Coefftrue            = [betatrue;SigsegTrue;rc_true];                % true parameters of the model
Coefficients         = zeros(nlin+nrcNL,1,4,R);
StandardErrors       = zeros(nlin+nrcNL,1,4,R);
ComputationTime      = zeros(1,1,4,R);

OwnElast             = zeros(nseg,nmkt,4,R);
CrossElastSameSeg    = zeros(nseg,nmkt,4,R);
CrossElastDiffSeg    = zeros(nseg,nmkt,4,R);
OneOwnElast          = zeros(1,nmkt,4,R);
OneCrossElastSameSeg = zeros(1,nmkt,4,R);
OneCrossElastDiffSeg = zeros(1,nmkt,4,R);

DiversionSameSeg    = zeros(nseg,nmkt,4,R);
DiversionDiffSeg    = zeros(nseg,nmkt,4,R);
OneDiversionSameSeg = zeros(1,nmkt,4,R);
OneDiversionDiffSeg = zeros(1,nmkt,4,R);

BiasOwnElast         = zeros(nseg,nmkt,3,R);
BiasCrossElastSameSeg= zeros(nseg,nmkt,3,R);
BiasCrossElastDiffSeg= zeros(nseg,nmkt,3,R);

FunctionValue        = zeros(1,1,4,R);
ExitFlag             = zeros(1,1,4,R);
AIC                  = zeros(1,1,4,R);
BIC                  = zeros(1,1,4,R);
% For the outside good I save the minimum, the average and the maximum for each simulation
OutsideShare         = zeros(3,R);
SaveSegment          = zeros(nobs,R);
SaveX1               = zeros(nobs,R);

for MC=1:R

    rng(seedvec(MC))
    
    % 3. Simulate data
    disp(['Simulation Numb: ',num2str(MC)])
    % Select the Seed
    
    % Generate X1 and Dstar according to a multivariate normal
    X1Dstar      = mvnrnd(meanc,Chelp1,nobs);
    
    Dstar        = X1Dstar(:,1);
    X1           = exp(X1Dstar(:,2));
    % The lognormal produces very extreme values sometimes. This
    % truncates but the distribution is not changed
    X1           = X1.*(X1 < 4) + (randn(nobs,1) + 3) .* (X1 > 4);
    multiX1      = reshape(X1,[prods,1,nmkt]);
    
    segment      = Dstar>cutoff;
    multisegment = reshape(segment,[prods,1,nmkt]);
    dumseg       = [Dstar<cutoff,Dstar>cutoff];
    multidumseg  = reshape(dumseg',[nseg,prods,nmkt]);
    multidumseg  = permute(multidumseg,[2 1 3]);
    multidumseg  = repmat(multidumseg,[1 1 1 Nnodes]);
    
    SaveSegment(:,MC) = segment;
    SaveX1(:,MC)      = X1;
    
    % random coefficient on X1
    xvu         = X1*nodedraws';
    multixvuL   = reshape(xvu,[prods,1,nmkt,Nnodes]);
    multixvu    = bsxfun(@times,multixvuL,multidumseg);
    
    muTrue      = xvu*rc_true;
    
    Xexo        = [const segment X1];
    deltaTrue   = Xexo*betatrue+ksi;
    
    % 4. Data Structure
    
    % Calculate the True/Observed Market shares
    Data.nmkt           = nmkt;
    Data.prods          = prods;
    Data.nobs           = nobs;
    Data.SigsegTrue     = SigsegTrue;
    Data.nseg           = nseg;
    Data.Nnodes         = Nnodes;
    Data.muTrue         = muTrue;
    Data.dumseg         = dumseg;
    Data.qweightrprods  = qweightrprods;
    Data.xvu            = xvu;
    Data.multidumseg    = multidumseg;
    Data.multixvu       = multixvu;
    Data.multixvuL      = multixvuL;
    Data.nlin           = nlin;
    Data.nrc            = nrc;
    Data.nrcNL          = nrcNL;
    Data.multiX1        = multiX1;
    Data.multisegment   = multisegment;
    
    % True Market Shares
    [~, ~, sj, sjg, s0, ~, ~, ~] = TrueShare(deltaTrue,Data);
    
    % 2D reshape
    sj                         = reshape(sum(sj,2),[prods*nmkt 1]);
    sjg                        = reshape(sum(sjg,2),[prods*nmkt 1]);
    s0                         = reshape(s0,[prods*nmkt 1]);
    y                          = log(sj) - log(s0);             % log-odds ratio
    lnsjg                      = log(sjg);
    lnsj                       = log(sj);
    
    OutsideShare(:,MC)         = [min(s0);mean(s0);max(s0)]; % save outside good shares (min, mean, max)
    
    Data.logobsshare           = lnsj;
    
    % 5 Optimal instruments
    
    % Compute optimal instruments
    % Delta without ksi
    PseudoDelta             = Xexo*betatrue;
    chamberlin              = jacobNL(NonLinTheta,PseudoDelta,Data);
    
    Z                       = [Xexo chamberlin];
    
    % 6 Estimation
    
    %% Logit
    % Weighting Matrix
    W           = (Z'*Z)\eye(size(Z,2));
    BLog        = ((Xexo'*Z*W*Z'*Xexo)\eye(size(Xexo,2)))*(Xexo'*Z*W*Z'*y);
    % Standard errors
    est         = y-Xexo*BLog;
    dgf         = (size(Xexo,1)-size(Xexo,2));
    ser         = (est'*est)./dgf;
    sst         = inv(Xexo'*Xexo);
    sterrLogit  = sqrt(ser*diag(sst));
    
    deltL       = y;
    alpha       = BLog(nlin);
    theta       = 0;
    % Logit elasticity wrt X1
    [ElaLogit, DiversionLogit]    = ElastLogit(alpha,theta,deltL,multiX1,Data);
    
    [MeanOwnElastLog,OneOwnElastLog,MeanCrossElastSameSegLog,OneCrossElastSameSegLog,MeanCrossElastDiffSegLog,OneCrossElastDiffSegLog,...
     MeanDivLog,MeanDivDiffLog,OneDivLog,OneDivDiffLog] = sumElast(ElaLogit,DiversionLogit,Data);
    
    % Function value
    FctValLog   = (est' * Z) * W * (Z' * est);
    % (size(Z),2) = number of moment conditions
    % size(BLog,1)= number of parameters
    % AIC = ln(nobs)*FunctValue-2*(numb moment conditions - numb parameters)
    AICLogit    = nobs * FctValLog - 2 * (size(Z,2)-size(BLog,1));
    % Add the two coefficients that are not estimated in the Logit (Nesting
    % parameter; Random Coefficient)
    % BIC = ln(nobs)*FunctValue-(numb moment conditions - numb parameters)*log(number of observations)
    BICLogit    = nobs * FctValLog - (size(Z,2)-size(BLog,1)) * log(nobs);
    
    % Store the results
    Coefficients(:,:,1,MC)          = [BLog;nan;nan];
    StandardErrors(:,:,1,MC)        = [sterrLogit;nan;nan];
    ComputationTime(:,:,1,MC)       = nan;
    OwnElast(:,:,1,MC)              = MeanOwnElastLog;
    
    CrossElastSameSeg(:,:,1,MC)     = MeanCrossElastSameSegLog;
 
    CrossElastDiffSeg(:,:,1,MC)     = MeanCrossElastDiffSegLog;

    OneOwnElast(:,:,1,MC)           = OneOwnElastLog;
    OneCrossElastSameSeg(:,:,1,MC)  = OneCrossElastSameSegLog;
    OneCrossElastDiffSeg(:,:,1,MC)  = OneCrossElastDiffSegLog;
       
    DiversionSameSeg(:,:,1,MC)     =  MeanDivLog;
    DiversionDiffSeg(:,:,1,MC)     =  MeanDivDiffLog;

    OneDiversionSameSeg(:,:,1,MC)     =  OneDivLog;
    OneDiversionDiffSeg(:,:,1,MC)     =  OneDivDiffLog;        
        
    FunctionValue(:,:,1,MC)         = FctValLog;
    ExitFlag(:,:,1,MC)              = nan;
    AIC(:,:,1,MC)                   = AICLogit;
    BIC(:,:,1,MC)               = BICLogit;
    
    
    %% Nested Logit
    XNestedLog  = [Xexo lnsjg];
    BNLog       = ((XNestedLog'*Z*W*Z'*XNestedLog)\eye(size(XNestedLog,2)))*(XNestedLog'*Z*W*Z'*y);
    % Standard errors
    est         = y-XNestedLog*BNLog;
    dgf         = (size(XNestedLog,1)-size(Xexo,2));
    ser         = (est'*est)./dgf;
    sst         = inv(XNestedLog'*XNestedLog);
    sterrNL     = sqrt(ser*diag(sst));
    deltNL     = y - BNLog(nlin+1)*lnsjg;
    % Nested Logit elasticity wrt X1
    alpha      = BNLog(nlin);
    theta      = [BNLog(nlin+1);0];
    [ElaNL, DiversionNL]     = ElastNestedLogit(alpha,theta,deltNL,multiX1,Data);
    
    [MeanOwnElastNL,OneOwnElastNL,MeanCrossElastSameSegNL,OneCrossElastSameSegNL,MeanCrossElastDiffSegNL,OneCrossElastDiffSegNL,...
     MeanDivNL,MeanDivDiffNL,OneDivNL,OneDivDiffNL] = sumElast(ElaNL,DiversionNL,Data);
    
    % Function value
    FctValNL   = (est' * Z) * W * (Z' * est);
    AICNLog    = nobs * FctValNL - 2 * (size(Z,2)-size(BNLog,1));
    BICNLog    = nobs * FctValNL - (size(Z,2)-size(BNLog,1)) * log(nobs);
    
    % Store the results
    Coefficients(:,:,2,MC)          = [BNLog;nan];
    StandardErrors(:,:,2,MC)        = [sterrNL;nan];
    ComputationTime(:,:,2,MC)       = nan;
    OwnElast(:,:,2,MC)              = MeanOwnElastNL;
    CrossElastSameSeg(:,:,2,MC)     = MeanCrossElastSameSegNL;
    CrossElastDiffSeg(:,:,2,MC)     = MeanCrossElastDiffSegNL;
    
    OneOwnElast(:,:,2,MC)           = OneOwnElastNL;
    OneCrossElastSameSeg(:,:,2,MC)  = OneCrossElastSameSegNL;
    OneCrossElastDiffSeg(:,:,2,MC)  = OneCrossElastDiffSegNL;
          
    DiversionSameSeg(:,:,2,MC)     =  MeanDivNL;
    DiversionDiffSeg(:,:,2,MC)     =  MeanDivDiffNL;

    OneDiversionSameSeg(:,:,2,MC)     =  OneDivNL;
    OneDiversionDiffSeg(:,:,2,MC)     =  OneDivDiffNL;
        
    FunctionValue(:,:,2,MC)         = FctValNL;
    ExitFlag(:,:,2,MC)              = nan;
    AIC(:,:,2,MC)                   = AICNLog;
    BIC(:,:,2,MC)                   = BICNLog;
    
    %% Random Coefficient Logit
    % Some data that speeds up computations
    xzwz            = Xexo'*Z*W*Z';
    Data.xzwz       = xzwz;
    xzwzx           = xzwz * Xexo;
    invxzwzx        = xzwzx\eye(size(xzwzx,2));
    Data.invxzwzx   = inv(xzwzx);
    Data.Z          = Z;
    Data.W          = W;
    Data.Xexo       = Xexo;
    % OPTIMIZATION
    % Options:see Knitro Option File KnitroBLP.opt
    % Pass Data into gmm using anonymous functions
    theta20RCL = randn(1,1);     % starting value
    delta0     = deltL;
    save mvalold delta0
    
    angmm    = @(theta20RCL)gmmLogit(theta20RCL,Data);
    
    options  = optimset( 'Display','iter',...
        'GradObj','on','TolCon',1E-6,...
        'TolFun',1E-6,'TolX',1E-6,...
        'Hessian', 'off','DerivativeCheck','off','FinDiffType','central');
    
    t1      = cputime;
    [thetaRCL, FctValRCL, exitflagRCL, ~, ~] = ...
        ktrlink(angmm,theta20RCL, [], [], [], [], [], [], [], options, 'knitroBLP.opt');
    
    % Ricalculate beta according to the thetaRCL that minimizes the function
    % given the use of multistart by KNITRO
    
    deltRCL      = deltaLogit(thetaRCL,Data);
    beta         = invxzwzx*(xzwz*deltRCL);
    th12RCL      = [beta ; thetaRCL];
    sterrRCL     = seRCLogit(th12RCL,Data);
    cputimegmm   = cputime-t1;
    
    % Random Coefficient Logit elasticity wrt X1
    alpha        = beta(nlin);
    [ElaRCLogit,DiversionRCLogit] = ElastLogit(alpha,thetaRCL,deltRCL,multiX1,Data);
    
    [MeanOwnElastRCL,OneOwnElastRCL,MeanCrossElastSameSegRCL,OneCrossElastSameSegRCL,MeanCrossElastDiffSegRCL,OneCrossElastDiffSegRCL,...
     MeanDivRCL,MeanDivDiffRCL,OneDivRCL,OneDivDiffRCL] = sumElast(ElaRCLogit,DiversionRCLogit,Data);
 
    AICRCLog     = nobs * FctValRCL - 2 * (size(Z,2)-size(th12RCL,1));
    BICRCLog     = nobs * FctValRCL - (size(Z,2)-size(th12RCL,1)) * log(nobs);
    
    
    % Store the results
    Coefficients(:,:,3,MC)          = [th12RCL(1:nlin);nan;abs(thetaRCL)];
    StandardErrors(:,:,3,MC)        = [sterrRCL(1:nlin);nan;sterrRCL(nlin+1)];
    ComputationTime(:,:,3,MC)       = cputimegmm;
    OwnElast(:,:,3,MC)              = MeanOwnElastRCL;
    CrossElastSameSeg(:,:,3,MC)     = MeanCrossElastSameSegRCL;
    CrossElastDiffSeg(:,:,3,MC)     = MeanCrossElastDiffSegRCL;
       
    OneOwnElast(:,:,3,MC)           = OneOwnElastRCL;
    OneCrossElastSameSeg(:,:,3,MC)  = OneCrossElastSameSegRCL;
    OneCrossElastDiffSeg(:,:,3,MC)  = OneCrossElastDiffSegRCL;
             
    DiversionSameSeg(:,:,3,MC)     =  MeanDivRCL;
    DiversionDiffSeg(:,:,3,MC)     =  MeanDivDiffRCL;

    OneDiversionSameSeg(:,:,3,MC)     =  OneDivRCL;
    OneDiversionDiffSeg(:,:,3,MC)     =  OneDivDiffRCL;
        
    FunctionValue(:,:,3,MC)         = FctValRCL;
    ExitFlag(:,:,3,MC)              = exitflagRCL;
    AIC(:,:,3,MC)                   = AICRCLog;
    BIC(:,:,3,MC)                   = BICRCLog;
    
    %% Random Coefficient Nested Logit
    
    t1      = cputime;
    
    % OPTIMIZATION
    % Options:see Knitro Option File KnitroBLP.opt
    % Pass Data into gmm using anonymous functions
    theta20RCNL  = [0.5;0.5];     % starting values - just pick a low value of pho to avoid overflow problems from the start
    delta0NL     = deltNL;
    save mvaloldNL delta0NL
    
    angmmNL = @(theta20RCNL)gmmNL(theta20RCNL,Data);
    
    options = optimset( 'Display','iter',...
        'GradObj','on','TolCon',1E-6,...
        'TolFun',1E-6,'TolX',1E-6,...
        'Hessian', 'off','DerivativeCheck','off','FinDiffType','central');
    
    %     [thetaRCNL, FctValRCNL, exitflagRCNL, output, gradient] = ...
    %         fminunc(angmmNL,theta20RCNL,options);
    
    x_L     = [0;-100];           % lower bound (only for the nested logit, to avoid numerical problems when pho gets negative)
    x_U     = [0.975;+100];       % upper bound (only for the nested logit, to avoid numerical problems when pho gets close to 1)
    [thetaRCNL, FctValRCNL, exitflagRCNL,~,~] = ...
        ktrlink(angmmNL,theta20RCNL, [], [], [], [], x_L, x_U, [], options, 'knitroBLP.opt');
    
    deltRCNL      = deltaNL(thetaRCNL,Data);
    betaNL        = invxzwzx*(xzwz*deltRCNL);
    th12RCNL      = [betaNL ; abs(thetaRCNL)];
    sterrRCNL     = seNL(th12RCNL,Data);
    cputimegmm    = cputime-t1;
    
    % Random Coefficient Nested Logit wrt X1
    alpha         = betaNL(nlin);
    
    [ElaRCNL,DiversionRCNL] = ElastNestedLogit(alpha,thetaRCNL,deltRCNL,multiX1,Data);
    
    [MeanOwnElastRCNL,OneOwnElastRCNL,MeanCrossElastSameSegRCNL,OneCrossElastSameSegRCNL,MeanCrossElastDiffSegRCNL,OneCrossElastDiffSegRCNL,...
     MeanDivRCNL,MeanDivDiffRCNL,OneDivRCNL,OneDivDiffRCNL] = sumElast(ElaRCNL,DiversionRCNL,Data);
 
    AICRCNLog    = nobs * FctValRCNL - 2 * (size(Z,2)-size(th12RCNL,1));
    BICRCNLog    = nobs * FctValRCNL - (size(Z,2)-size(th12RCNL,1)) * log(nobs);
    
    % Store the results
    Coefficients(:,:,4,MC)          = th12RCNL;
    StandardErrors(:,:,4,MC)        = sterrRCNL;
    ComputationTime(:,:,4,MC)       = cputimegmm;
    OwnElast(:,:,4,MC)              = MeanOwnElastRCNL;
    CrossElastSameSeg(:,:,4,MC)     = MeanCrossElastSameSegRCNL;
    CrossElastDiffSeg(:,:,4,MC)     = MeanCrossElastDiffSegRCNL;
               
    OneOwnElast(:,:,4,MC)           = OneOwnElastRCNL;
    OneCrossElastSameSeg(:,:,4,MC)  = OneCrossElastSameSegRCNL;
    OneCrossElastDiffSeg(:,:,4,MC)  = OneCrossElastDiffSegRCNL;
    
    DiversionSameSeg(:,:,4,MC)     =  MeanDivRCNL;
    DiversionDiffSeg(:,:,4,MC)     =  MeanDivDiffRCNL;

    OneDiversionSameSeg(:,:,4,MC)     =  OneDivRCNL;
    OneDiversionDiffSeg(:,:,4,MC)     =  OneDivDiffRCNL;
        
    FunctionValue(:,:,4,MC)         = FctValRCNL;
    ExitFlag(:,:,4,MC)              = exitflagRCNL;
    AIC(:,:,4,MC)                   = AICRCNLog;
    BIC(:,:,4,MC)                   = BICRCNLog;
    
    savefile = '1resultsRCNL_LowCorr_Pho03Sigma1_Dstar0.mat';
    
    save    (savefile, 'Coefftrue', 'Coefficients', 'StandardErrors', 'ComputationTime',...
        'OwnElast', 'CrossElastSameSeg', 'CrossElastDiffSeg',...
        'OneOwnElast','OneCrossElastSameSeg','OneCrossElastDiffSeg',...
        'DiversionSameSeg', 'DiversionDiffSeg', 'OneDiversionSameSeg','OneDiversionDiffSeg',...
        'BiasOwnElast', 'BiasCrossElastSameSeg', 'BiasCrossElastDiffSeg',...
        'FunctionValue', 'ExitFlag', 'AIC', 'BIC', ...
        'SaveSegment', 'SaveX1', 'OutsideShare', 'corrDstarX1')
    
end