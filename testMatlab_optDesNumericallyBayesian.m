%% Specify individual or distribution tissue propeties
tissueProperties='distribution';%'distribution';%'individual'

%% Specify clinical or preclinical protocol
clin_or_preclin='preclinical';%'clinical' or 'preclinical'

%% Constrain TE and bmax (only for preclinical for now)
te_bmax_constrain='constrainTEbmax';

%% Generate tissue parameter matrix -
%  samples from 4D prior distribution for robust design
switch tissueProperties
    case 'individual'
        radVal=10e-6;%/(1/0.4^(1/3));
        diVal=1.5e-9;%1e-9;
        deVal=2e-9;%2e-9;1.5e-9
        fiVal=0.60;
        t2Val=90e-3;
        
        tissueParamMatrix=[radVal diVal deVal fiVal t2Val];
        
    case 'distribution'
        
        numSamples=100;%50;%
        
        %t2=90e-3; %~90 ms: prostate tumour median T2, Langer et al. JMRI 2009
        t2=65e-3; % For phantom optimised protocol
        
        radMin=5e-6;
        radMax=20e-6;
        
        diMin=0.5e-9;
        diMax=3e-9;
        
        deMin=0.5e-9;
        deMax=3e-9;
        
        fiMin=0.1;
        fiMax=0.74;
        
        RandStream.setGlobalStream...
            (RandStream('mt19937ar','seed',1748));
        tissueParamMatrix=zeros(numSamples,4);
        for i=1:numSamples
            tissueParamMatrix(i,:)=...
                unifrnd([radMin diMin deMin fiMin],[radMax diMax deMax fiMax]);
        end
        tissueParamMatrix=cat(2,tissueParamMatrix,repmat(t2,numSamples,1));
        
end
%% SNR options
snr=29;%50;%20;%
snrLowLim=3;

%% Seeds
seedList=[231187 910];%[11160 231187 910];

%% GA population size
populationSize=1600; %200,400,800,1600

%% Derivative spacing factor - if directory name doesn't include this, 
%  a factor of 5 will have been used
%  01/12/16 - value was included in directory name when testing a value of 
%  8, but now carrying on with 5, so remove this suffix for consistency 
%  with previous names 
derivSpacingFactor=5;

%% Make directory
switch tissueProperties
    case 'individual'
        newDir=strcat('scaled_bayesianDOptGENETICALGORITHM_',...
            clin_or_preclin,'_',...
            te_bmax_constrain,'_',...
            'rVal_',num2str(radVal*1e6),'_',...
            'diVal_',num2str(diVal*1e9),'_',...
            'deVal_',num2str(deVal*1e9),'_',...
            'fiVal_',num2str(fiVal*100),'_',...
            'T2Val_',num2str(t2Val*1e3),'_',...
            'snr_B0',num2str(snr),'_',...
            'snrLowerLim',num2str(snrLowLim),'_',...       
            'gaPopulation',num2str(populationSize));
            %'derivSpacingFactor',num2str(derivSpacingFactor));
        if ~exist(newDir, 'dir')
            mkdir(newDir)
        else
            error('!!! Directory for these tissue property ranges already exists')
        end
        
    case 'distribution'
        
        newDir=strcat('scaled_bayesianDOptGENETICALGORITHM_',...
            clin_or_preclin,'_',...
            te_bmax_constrain,'_',...
            'rRange_',num2str(radMin*1e6),'_',num2str(radMax*1e6),'_',...
            'diRange_',num2str(diMin*1e9),'_',num2str(diMax*1e9),'_',...
            'deRange_',num2str(deMin*1e9),'_',num2str(deMax*1e9),'_',...
            'fiRange_',num2str(fiMin*100),'_',num2str(fiMax*100),'_',...
            'T2',num2str(t2*1e3),'_',...
            'snr_B0',num2str(snr),'_',...
            'snrLowerLim',num2str(snrLowLim),'_',...
            'numSamples',num2str(numSamples),'_',...
            'gaPopulation',num2str(populationSize));
            %'derivSpacingFactor',num2str(derivSpacingFactor));
        if ~exist(newDir, 'dir')
            mkdir(newDir)
        else
            error('!!! Directory for these tissue property ranges already exists')
        end
        
end

%% Call optimisation function and save output
%{%
%% Loop over seeds
for j=1:numel(seedList)
    disp(j)
    structfield=strcat('seed',num2str(seedList(j)));
    optDesStruct.(structfield)=...
        optDesNumericalObjFnBayesian_DiDeT2('murday_cotts',...
        tissueParamMatrix,snr,snrLowLim,...
        [0 0 0 0; 0 0 0 0; 0 0 0 0],...
        [Inf Inf Inf Inf ;Inf Inf Inf Inf; Inf Inf Inf Inf],...
        seedList(j),populationSize,derivSpacingFactor,clin_or_preclin,...
        te_bmax_constrain);
    try
        optDesStruct.(structfield).Dopt
    catch
    end
    % Save output as we go, in case matlab crashes
    save(strcat(newDir,'/optDesStruct.mat'),'optDesStruct')
end
%}

%% Save nonOpt design and tissueProp vector, for 'individual' directories
%  These are needed for fitting simulations (see
%  ./fittingSimulations/fittingScript_cf_nonOpt_localOpt_robustOpt.m)
switch tissueProperties
    case 'individual'
        switch clin_or_preclin
            case 'clinical'
                gradNonOpt=[0.03;0.06]';
                DELNonOpt=[0.02;0.072]';
                delNonOpt=0.015;
            case 'preclinical'
                gradNonOpt=[0.093;0.187;0.280]';
                DELNonOpt=[0.012;0.025;0.045]';
                delNonOpt=0.004;
        end
        design=cat(2,combvec(gradNonOpt,DELNonOpt)',...
            repmat(delNonOpt,size(combvec(gradNonOpt,DELNonOpt)',1),1));
        design=design.*1e3; % Designs are scaled now!
        save(strcat(newDir,'/nonOpt.mat'),'design')
        
        save(strcat(newDir,'/tissueProp.mat'),'tissueParamMatrix')
        
    case 'distribution'
        disp('For distribution directories, nonOpt and tissueProp are not needed')
end