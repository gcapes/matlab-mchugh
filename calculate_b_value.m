function b = calculate_b_value(N, G, d, s, D)

% N = no of oscillations; N=0 standard PGSE; gr_repetitions/2
% G = max gradient amplitude in mT/m; gr_str_factor_max
% d = total duration of each osc gradient in ms; gr_dur * gr_repetitions
% s = duration of slope of gradient in ms; gr_slope
% D = big delta = separation of gradients (only for N=0); 

gamma = 2*pi*42.57746778e6;%267.5e6; % gyromagnetic ratio in rad s-1 T-1
G = G/1e6; %in T/mm
s = s/1000; % in s

if N==0
    D = D/1000; %in s
    delta = d/1000 - s; % in s
    b  = gamma.^2 .* G.^2 .* (delta.^2 .* (D - delta./3) + s.^3./30 - delta.*(s.^2)./6);
    
else

    delta = d./(2.*1000.*N) - s;  %in s
    b = gamma.^2 .* G.^2 .* N .*(4.*delta.^3./3 + 2.*s.*delta.^2 - delta.*s.^2./3 + s.^3./15);
end

%G_achieved=sqrt(round(b)./(gamma.^2 .* (delta.^2 .* (D - delta./3) + s.^3./30 - delta.*(s.^2)./6))).*1e6; %in mT/m
%((G.*1e6-G_achieved)./(G.*1e6)).*100
%b-value for fourier friendly waveform
% k=d/(4*1000*N) - s;
% b_ff=2*(gamma^2)*G^2*( (2*N-1)*( delta^3/3 + s*delta^2/2 - delta*s^2/12 + s^3/60) + 2*(k^3/3 +s*k^2/2 - k*s^2/12 + s^3/60) );
% b_ff_OVER_b=b_ff/b;
% plateau=d/(4*1000*N)- 2*s;

%D_eff = D - (delta/3);
%k = 4*pi*pi*D_eff;
%q_sq = b./k;
%q = sqrt(q_sq)*10; %gives q value in cm^-1