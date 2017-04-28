function [optDesNumericalObjFnBayesian_DiDeT2] = ...
    optDesNumericalObjFnBayesian_DiDeT2(sig_model,tissueParamSamples,snrB0,...
    snrLowerLim,lwrBnds,upprBnds,initialParamsSeed,populationSize,...
    derivativeSpacingFactor,clin_or_preclin,constrain_TE_bmax)

%% Inputs:
%  sig_model - signal model name (currently only works for 'murday_cotts')
%  tissueParamSamples - nx4 matrix, with n values of 4 tissue properties
%                       (radius, intracellular diffusivity,
%                       extracellular diffusivity, volume fraction)
%                     - acts as prior distribution for Bayesian/robust
%                       design
%  snrB0 - signal to noise ratio (in b=0 data)
%  snrLowerLim - lower SNR limit, to avoid bias due to Rician noise
%  lwrBnds - lower bounds on optimised G, DELTA, delta values
%  upprBnds - upper bounds on optimised G, DELTA, delta values
%  initialParamsSeed - seed for RNG
%  populationSize - population size for GA
%  derivativeSpacingFactor - order of magnitude for step size in derivative
%                            calculation
%  clin_or_preclin - 'clinical' or 'preclinical' specifying scanner
%  constraints
%  constrain_TE_bmax - 'constrainTEbmax' to constrain TE and bmax to hard coded values

%% Output:
%  optDesNumericalObjFnBayesian_DiDeT2 - structure containing results and
%                                        information about minimisation
%                                      - Dopt field contains optimum G,
%                                        DELTA and delta values (4x3 matrix)

%% Same as optDesNumericalObjFn_DiDeT2, but calculates one Bayesian/robust
%  design, based on sampling a prior distribution of tissue parameters
%  - see p. 306 in Atkinson 'Optimium Experimental Designs, with SAS'.

%% Do we need to scale parameters or not?
%  05/05/16 - SCALING APPEARS TO HELP - WARNINGS ARE NOT GIVEN WHEN USING
%  SCALED PARAMETERS. PRESUMABLY RELATED TO MATRIX INVERSION & DETERMINANT?
scaling='scaled'; %% 'scaled' or 'unscaled'

switch scaling
    case 'scaled' % if scaling, make sure tissueParams are scaled
        tissueParamSamples=tissueParamSamples.*repmat([1e6 1e9 1e9 1 1e3],...
            size(tissueParamSamples,1),1);
        size(tissueParamSamples)
    case 'unscaled'
        % no need to do anything to the parameters
end

%% Noise standard deviation from snrB0
%  snrB0 assumed to be for TE = 21.3 ms and T2 given in tissueParamSamples
snr_ref=snrB0;
TE_ref=21.3e-3;
switch scaling
    case 'scaled'
        TE_ref=TE_ref*1e3;
    case 'unscaled'
        % no need to do anything 
end
noiseSD=exp(-TE_ref/tissueParamSamples(1,5))/snr_ref;

%% Choose tissue parameter combination to use for SNR constraint calculations
%  Take median for all parameters for now - may need re-thinking!
tissueParamUseForConstraints=median(tissueParamSamples,1);

%% Set RNG
RandStream.setGlobalStream...
    (RandStream('mt19937ar','seed',initialParamsSeed));

%% Load roots for MC equation
rootsS=loadRename('uniqueRootsS_first65_higherPrecision.mat');

%% Generate starting values, depending on signal model
switch sig_model
    case 'murday_cotts'
        % Starting value matrix
        numStartingVals=populationSize;
        
        switch clin_or_preclin
            case 'clinical'
                te=Inf; %seconds
                readoutAndPreFirstGradTime=13e-3; %seconds; see timings in books (22/12/14 & 5/9/12)
                %timeFor180AndCrushers=5e-3; %see pulse sequence diagram book 4, 13/9/11
                timeFor180AndCrushers=12e-3; %180816 - previously 5ms, but needs increasing to 12ms for current 3T protocols
                gmax=60e-3; %T/m
                bValMin=150; %s/mm^2
                % Make bvalMax very high but not Inf so @nonlin_constraints doesn't
                % return Inf
                bValMax=1e9; %s/mm^2
                %bValMax=Inf; %s/mm^2
                protonGamma=2*pi*42.57746778e6; %proton gyromagnetic ratio
                gList=(0:5:65).*1e-3;
                delList=(1:5:100).*1e-3;
                DELList=(1:5:200).*1e-3;
            case 'preclinical'
                te=Inf; %seconds
                readoutAndPreFirstGradTime=13e-3; %seconds; see timings in books (22/12/14 & 5/9/12)
                timeFor180AndCrushers=5e-3; %seconds; see pulse sequence diagram book 4, 13/9/11
                gmax=300e-3; %T/m
                bValMin=150; %s/mm^2
                bValMax=1e9; %s/mm^2
                protonGamma=2*pi*42.57746778e6; %proton gyromagnetic ratio
                gList=(0:10:300).*1e-3;
                delList=(1:10:100).*1e-3;
                DELList=(1:10:200).*1e-3;
                switch constrain_TE_bmax
                    case 'constrainTEbmax'
                        disp('!!! Constraining TE and bmax !!!')
                        te=55e-3;
                        bValMax=1950;
                    otherwise
                        % Leave values as above
                end     
        end
        switch scaling
            case 'unscaled'
                % no need to do anything to the parameters above
            case 'scaled'
                te=te.*1e3; %milliseconds
                readoutAndPreFirstGradTime=readoutAndPreFirstGradTime.*1e3; %milliseconds; see timings in books (22/12/14 & 5/9/12)
                timeFor180AndCrushers=timeFor180AndCrushers.*1e3; %milliseconds; see pulse sequence diagram book 4, 13/9/11
                gmax=(gmax./1e9).*1e12;
                % Be careful with b-value scaling, as for this calculation
                % values get scaled in generateScheme calculations and do
                % not scale protonGamma, whereas for signal calculation, 
                % protonGamma IS scaled and scale G,DEL,del are used
                % directly
                bValMin=(bValMin*1e15);%(bValMin/1e9); %
                bValMax=(bValMax*1e15);%(bValMax/1e9); %
                protonGamma=2*pi*42.57746778e6; %proton gyromagnetic ratio
                protonGamma=protonGamma./1e12; % needs scaling as it is passed to MCanalytic_Di_De_WITH_T2_passRoots
                gList=(gList./1e9).*1e12;
                delList=delList.*1e3;
                DELList=DELList.*1e3;
        end
        
        %% Generate scan parameters and remove those which result in snr<snrLowerLim
        scanParamsAll=generateSchemes(gList,DELList,delList,protonGamma,...
            timeFor180AndCrushers,gmax,te,readoutAndPreFirstGradTime...
            ,bValMin,bValMax,1);
        minTE=scanParamsAll(:,2)+scanParamsAll(:,3)+...
            readoutAndPreFirstGradTime;
        scanParamsAll=cat(2,scanParamsAll,minTE);   
        allSigs=MCanalytic_Di_De_WITH_T2_passRoots(...
            tissueParamUseForConstraints,scanParamsAll,rootsS);
        ind=allSigs/noiseSD>snrLowerLim;
        scanParams=scanParamsAll(ind,:);
                
        %% Select subset of possible combinations
        for strtVal=1:numStartingVals
            randInd=randi(size(scanParams,1),4,1);
            startingVals(:,:,strtVal)=scanParams(randInd,1:3);
        end
        
    otherwise
        error('!!! invalid sig_model')
end

%% Specify linear constraints (non-linear constraints dealt with below)
%lincon1=fitParams(3)+timeFor180AndCrushers-fitParams(2);
%lincon2=fitParams(3)+fitParams(2)-te+readoutAndPreFirstGradTime;
%lincon3=fitParams(1)-gmax;
%nonlincon4=bValMin-calculate_b_value(0,fitParams(1).*1e3,fitParams(3).*1e3,0,fitParams(2).*1e3);
%nonlincon5=calculate_b_value(0,fitParams(1).*1e3,fitParams(3).*1e3,0,fitParams(2).*1e3)-bValMax;
%nonlincon6=snrLowerLim - (MCanalytic_Di_De_WITH_T2_passRoots(tissueParam,...
%    [fitParams(:,1) fitParams(:,2) fitParams(:,3) ...
%    repmat(2*pi*42.57746778e6,size(fitParams,1),1) ...
%            fitParams(:,2)+fitParams(:,3)+readoutAndPreFirstGradTime])/noiseSD)
linearConstraintMatrix=...
   [0 0 0 0 -1 0 0 0 1 0 0 0;
    0 0 0 0 0 -1 0 0 0 1 0 0;
    0 0 0 0 0 0 -1 0 0 0 1 0;
    0 0 0 0 0 0 0 -1 0 0 0 1;
    0 0 0 0 1 0 0 0 1 0 0 0;
    0 0 0 0 0 1 0 0 0 1 0 0;
    0 0 0 0 0 0 1 0 0 0 1 0;
    0 0 0 0 0 0 0 1 0 0 0 1;
    1 0 0 0 0 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0 0 0 0];
linearConstraintVec=...
    [repmat(-timeFor180AndCrushers,size(startingVals,1),1);...
    repmat(te-readoutAndPreFirstGradTime,size(startingVals,1),1);...
    repmat(gmax,size(startingVals,1),1)];

switch sig_model
    case 'murday_cotts'
        %{
            [fitParams,objFnVal,~,~]=...
                fminsearchcon(@optDesNumericalObFnBayesian,...
                startingVals(:,:,strtValInd),lwrBnds,upprBnds,...
                linearConstraintMatrix,linearConstraintVec,...
                @nonlin_constraints,options);
            ofv(strtValInd+1)=objFnVal;
            % Update estimates if obj fun is lower
            if ofv(strtValInd+1)<min(ofv(1:strtValInd))
                Gfit=fitParams(:,1);
                DELfit=fitParams(:,2);
                delfit=fitParams(:,3);
                strtValIndAccepted=strtValInd;
                finalObjFnVal=ofv(strtValInd+1);
            else
            end
        %}
        %{
            [fitParams,objFnVal,extflg,output]=...
             fmincon(@optDesNumericalObFnBayesian,...
                startingVals(:,:,strtValInd),...
                linearConstraintMatrix,linearConstraintVector,...
                [],[],lwrBnds,upprBnds,...
                @nonlin_constraints,options);
        %}
        %{%
        %% Starting values - organise as vector: 
        %  G,G,G,G,DEL,DEL,DEL,DEL,del,del,del,del
        for ii=1:size(startingVals,3)
            thisMat=startingVals(:,:,ii);
            X0(ii,:)=thisMat(:)';
        end
       
        %% Run Genetic Algorithm
        options = gaoptimset('InitialPopulation',X0,...
            'MutationFcn',@mutationadaptfeasible,'Display','iter',...
            'StallGenLimit',100,'PopulationSize',populationSize);
            %'PlotFcns',@gaplotbestf,'Vectorized','on');
        [fitParams,objFnVal,extflg,output,population]=...
            ga(@optDesNumericalObFnBayesian,...
            12,...
            linearConstraintMatrix,linearConstraintVec,...
            [],[],lwrBnds,upprBnds,...
            @nonlin_constraints,options);
        %% Assign output
        fitParams=reshape(fitParams,...
            size(startingVals,1),size(startingVals,2));
        Gfit=fitParams(:,1);
        DELfit=fitParams(:,2);
        delfit=fitParams(:,3);
        finalObjFnVal=objFnVal;
end

%% See if obj fn can be improved by searching in region of minima found above
%% NOT USED FOR GENETIC ALGORITHM
%{
% Perturb values by a given percentage
pc_perturb=10;
startingVals_=startingVals(:,:,strtValIndAccepted);
startingVals_=startingVals_(:);

% Start at lowest obj fn value found above
ofvNew=finalObjFnVal;

% Assign potential new outputs as empty arrays
Gfit_second_round=[];
DELfit_second_round=[];
delfit_second_round=[];
strtValIndAccepted_second_round=[];
finalObjFnVal_second_round=ofvNew;
for ii=1:numel(startingVals_)*4
    perturbStartingVals=startingVals_;
    if ismember(ii,1:numel(startingVals_))
        disp('Starting value perturbation - adding 10%')
        perturbStartingVals(ii)=...
            perturbStartingVals(ii)+perturbStartingVals(ii)*(pc_perturb/100);
    elseif ismember(ii,numel(startingVals_)+1:2*numel(startingVals_))
        disp('Starting value perturbation - subtracting 10%')
        ind=ii-numel(startingVals_);
        perturbStartingVals(ind)=perturbStartingVals(ind)-...
            perturbStartingVals(ind)*(pc_perturb/100);
    elseif ismember(ii,2*numel(startingVals_)+1:3*numel(startingVals_))
        disp('Starting value perturbation - adding 5%')
        ind=ii-2*numel(startingVals_);
        perturbStartingVals(ind)=perturbStartingVals(ind)+...
            perturbStartingVals(ind)*(pc_perturb/200);
    elseif ismember(ii,3*numel(startingVals_)+1:4*numel(startingVals_))
        disp('Starting value perturbation - subtracting 5%')
        ind=ii-3*numel(startingVals_);
        perturbStartingVals(ind)=perturbStartingVals(ind)-...
            perturbStartingVals(ind)*(pc_perturb/200);
    end
    
    perturbStartingVals=reshape(perturbStartingVals,size(startingVals,1),...
        size(startingVals,2));
    
    % Check that perturbed values satisfy linear and nonlinear constraints
    % - if not, skip this perturbation
    nlinconCheck=nonlin_constraints(perturbStartingVals);
    linconCheck=linearConstraintMatrix*perturbStartingVals(:);
    if any(nlinconCheck>0) || any(linconCheck>linearConstraintVec)
        ofvNew(ii+1)=NaN;
        continue
    else
    end
    
%{
    switch scaling
        case 'scaled'
    [perturbStartingVals(:,1).*1e12 perturbStartingVals(:,2) perturbStartingVals(:,3)]
        case 'unscaled'
    [perturbStartingVals(:,1) perturbStartingVals(:,2) perturbStartingVals(:,3)]
    end
%}
    
    % Minimisation
    switch sig_model
        case 'murday_cotts'
            [fitParams,objFnVal,~,~]=...
                fminsearchcon(@optDesNumericalObFnBayesian,...
                perturbStartingVals,lwrBnds,upprBnds,...
                linearConstraintMatrix,linearConstraintVec,...
                @nonlin_constraints,options);
            
            ofvNew(ii+1)=objFnVal;
            % Update estimates if obj fun is lower
            if ofvNew(ii+1)<min(ofvNew(1:ii))
                Gfit_second_round=fitParams(:,1);
                DELfit_second_round=fitParams(:,2);
                delfit_second_round=fitParams(:,3);
                strtValIndAccepted_second_round=ii;
                finalObjFnVal_second_round=ofvNew(ii+1);
            else
            end
    end
end
%}

%% Collect output
switch sig_model
    case 'murday_cotts'
        
        fitParamsAccepted.scanParams=cat(2,Gfit(:),DELfit(:),delfit(:));
        fitParamsAccepted.finalObjFnVal=finalObjFnVal;
        fitParamsAccepted.exitFlag=extflg;
        fitParamsAccepted.population=population;
        fitParamsAccepted.Dopt=fitParamsAccepted.scanParams;
        
        [fitParamsAccepted.Dopt(:,1)...
            fitParamsAccepted.Dopt(:,2) ...
            fitParamsAccepted.Dopt(:,3)]    
        %% NOT USED FOR GENETIC ALGORITHM
        %{
        fitParamsAccepted.scanParams_second_round=...
            cat(2,Gfit_second_round(:),DELfit_second_round(:),delfit_second_round(:));
        fitParamsAccepted.strtValAccepted_second_round=strtValIndAccepted_second_round;
        fitParamsAccepted.finalObjFnVal_second_round=finalObjFnVal_second_round;
        
        if fitParamsAccepted.finalObjFnVal_second_round<...
                fitParamsAccepted.finalObjFnVal
            fitParamsAccepted.Dopt=...
                fitParamsAccepted.scanParams_second_round;
        else
            fitParamsAccepted.Dopt=fitParamsAccepted.scanParams;
        end
        %}
        optDesNumericalObjFnBayesian_DiDeT2=fitParamsAccepted;
        
end

%% Non-linear constraints
    function [c,ceq]=nonlin_constraints(fitParams)
        switch sig_model
            case 'murday_cotts'
                %% Re-organise vector to matrix
                fitParams=reshape(fitParams,...
                    size(startingVals,1),size(startingVals,2));
                c=[bValMin-calculate_b_value(0,fitParams(:,1).*1e3,fitParams(:,3).*1e3,0,fitParams(:,2).*1e3);...
                    calculate_b_value(0,fitParams(:,1).*1e3,fitParams(:,3).*1e3,0,fitParams(:,2).*1e3)-bValMax;...
                    snrLowerLim - (MCanalytic_Di_De_WITH_T2_passRoots(tissueParamUseForConstraints,...
                    [fitParams(:,1) fitParams(:,2) fitParams(:,3) ...
                    repmat(protonGamma,size(fitParams,1),1) ...
                    fitParams(:,2)+fitParams(:,3)+readoutAndPreFirstGradTime],rootsS)/noiseSD)
                    ];
                ceq=[];
            otherwise
                error('!!! invalid sig_model')
        end
        
    end

    function x=optDesNumericalObFnBayesian(a)
        switch sig_model
            case 'murday_cotts'
                %% Re-organise vector to matrix
                aTemp=[];
                for icnt=1:size(a,1)
                    thisVec=a(icnt,:);
                    aTemp=cat(1,aTemp,reshape(thisVec,...
                        size(startingVals,1),size(startingVals,2)));
                end
                a=aTemp;
                
                %% Loop over sampled tissue parameters
                neglogdetI=zeros(1,size(tissueParamSamples,1)); %preallocate
                for tissueSample=1:size(tissueParamSamples,1)
                    tissueParam=tissueParamSamples(tissueSample,:);
                    
                    %% First get derivatives of signal model for each tissueParam
                    spacingStartList=10.^(floor(log10(abs([...
                        tissueParam(1) tissueParam(2) tissueParam(3) tissueParam(4)])))); %param order of mag
                    spacingList=spacingStartList./(10^(derivativeSpacingFactor));
                    collectSigDeriv=...
                        zeros(size(a,1),numel(tissueParam)-1);
                    for i=1:numel(tissueParam)-1 %Just include R, Di, De, fi
                        modelParamMat=[];
                        
                        %get spacing for this model parameter
                        spacing=spacingList(i);
                        
                        %get model parameter to differentiate wrt, and vector of values to
                        %calculate numerical gradient with
                        modelParamDiff_wrt=tissueParam(i);
                        modelParamDiff_wrt_vec=...
                            ((modelParamDiff_wrt-spacing):spacing:(modelParamDiff_wrt+spacing))';
                        
                        %for the other parameters, set them to the values specified above
                        [~,modelParamsSet]=...
                            find(tissueParam.*ismember(1:numel(tissueParam),i)==0);
                        modelParamMat(modelParamsSet)=tissueParam(modelParamsSet);
                        modelParamMat=repmat(modelParamMat,numel(modelParamDiff_wrt_vec),1);
                        modelParamMat(:,i)=modelParamDiff_wrt_vec;
                        
                        %loop over modelParamMat
                        thisScanMat=[a(:,1) a(:,2) a(:,3) ...
                            repmat(protonGamma,size(a,1),1) ...
                            a(:,2)+a(:,3)+readoutAndPreFirstGradTime];
                        for j=1:size(thisScanMat,1)
                            sig=zeros(1,size(modelParamMat,1));
                            thisScan=thisScanMat(j,:);
                            for k=1:size(modelParamMat,1)
                                thisModel=modelParamMat(k,:);
                                sig(k)=MCanalytic_Di_De_WITH_T2_passRoots(...
                                    thisModel,thisScan,rootsS);
                            end
                            %
                            %for 1x3 vector (i.e. one value either side of specified value),
                            %gradient will be calculated using central difference (see code for
                            %built-in gradient function)
                            sigDeriv=gradient(sig,spacing);
                            
                            %find derivative corresponding to specified model parameter
                            collectSigDeriv(j,i)=...
                                sigDeriv(modelParamDiff_wrt_vec==modelParamDiff_wrt);
                        end
                    end
                    
                    % calculate information matrix and summary statistic
                    % for this tissueSample
                    F=collectSigDeriv;
                    I=F'*F;
                    neglogdetI(tissueSample)=-log(det(I));
                end
            otherwise
                error('!!! invalid sig_model')
        end
        
        %% Minimise sum over x values for all tissue samples
        x=sum(neglogdetI);
    end
end