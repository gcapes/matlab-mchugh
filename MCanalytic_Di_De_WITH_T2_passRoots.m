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
for m=1:numel(rootsS)
    aS=1/(alpS(m).^2.*(alpS(m).^2.*R.^2-2));
    bS=(2.*del)./(alpS(m).^2.*Di);
    cS=(2+exp(-alpS(m).^2.*Di.*(DEL-del)) - 2.*exp(-alpS(m).^2.*Di.*del) -2.*exp(-alpS(m).^2.*Di.*DEL) + exp(-alpS(m).^2.*Di.*(DEL+del)) )./((alpS(m).^2.*Di).^2);
    exprM=aS.*(bS-cS);
    sumM=sumM+exprM;
end
murdayCotts=exp(-2.*(gamma.^2).*(G.^2).*sumM);
F=S0.*exp(-TE./T2).*((f.*murdayCotts)+((1-f).*exp(-(((G.*del.*gamma).^2).*(DEL-del./3)).*(De./(1+f./2)))));
end