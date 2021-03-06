function F=MCanalytic_Di_De_WITH_T2_passRoots(a,b,rootsS)
%function to generate signals from impermeable model for specific tissue
%properties (a) and specific pulse sequence parameters (b)

% -identical to MCanalytic_Di_De_WITH_T2 but allows equation roots to be
% passed directly, so that the values don't have to be loaded on every call
% i.e. rootsS=loadRename('uniqueRootsS_first65_higherPrecision.mat');

sumM=0;
%%
%tissue properties (assume S0=1)
R=a(1);
Di=a(2);
De=a(3);
f=a(4);
T2=a(5);
S0=1;

%PGSE scan parameters
G=b(:,1);
DEL=b(:,2);
del=b(:,3);
gamma=b(:,4);
TE=b(:,5);

%calculate signal from two-compartment analytic expression
alpS=rootsS./R;

nroot=numel(rootsS);
ndel=numel(del);

% Vectorize assignments
aSvec = 1./(alpS.^2 .*(alpS.^2 .*R.^2-2));
bSvec = repmat(2*del,1,nroot)./repmat((alpS.^2.*Di)',ndel,1);

% Duplicate vectors into matrices to enable vectorised assignment
alpS_vsq=repmat(alpS,1,ndel)'.^2;
del_vec=repmat(del,1,nroot);
DEL_vec=repmat(DEL,1,nroot);
cSvec = (2+exp(-alpS_vsq.*Di.*(DEL_vec-del_vec)) - 2.*exp(-alpS_vsq.*Di.*del_vec) -2.*exp(-alpS_vsq.*Di.*DEL_vec) + exp(-alpS_vsq.*Di.*(DEL_vec+del_vec)) )./((alpS_vsq.*Di).^2);

aSvec_vec=repmat(aSvec,1,ndel)';
exprM_vec=aSvec_vec.*(bSvec-cSvec);

sumMvec = sum(exprM_vec,2);

murdayCotts=exp(-2.*(gamma.^2).*(G.^2).*sumMvec);
F=S0.*exp(-TE./T2).*((f.*murdayCotts)+((1-f).*exp(-(((G.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+f./2)))));
end
