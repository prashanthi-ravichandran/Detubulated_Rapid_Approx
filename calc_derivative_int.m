function [dCaiArray, dCaNSRArray,dLTRPNCaArray, dHTRPNCaArray,  JMyo_diff, JNSR_diff, Jup, Jtrpn] = calc_derivative_int(CaiArray, CaNSRArray, LTRPNCaArray, HTRPNCaArray, Jtr_int, JRyR_int, diff_factor)

global r Vmyo VNSR n_shells dr 

tol = 1e-6;

% Buffering parameters
% total troponin low affinity site conc. (mM)
LTRPNtot= 70.0e-3;
% total troponin high affinity site conc. (mM)
HTRPNtot= 140.0e-3;
% Ca++ on rate for troponin high affinity sites (1/(mM*ms))
khtrpn_plus= 20.0;
% Ca++ off rate for troponin high affinity sites (1/ms)
khtrpn_minus= 0.066e-3;
% Ca++ on rate for troponin low affinity sites (1/(mM*ms))
kltrpn_plus= 40.0;
% Ca++ off rate for troponin low affinity sites (1/ms)
kltrpn_minus= 40.0e-3;
% total myoplasmic calmodulin concentration (mM)
CMDNtot= 50.0e-3;
% total myoplasmic EGTA concentration (mM)
EGTAtot= 0.0;
% Ca++ half sat. constant for calmodulin (mM)
KmCMDN= 2.38e-3;
% Ca++ half sat. constant for EGTA (mM)
KmEGTA= 1.5e-4;

VCSR= 2.777042242067234e-11; % network SR volume (uL)
% SR parameters
Kfb=0.26e-3; % foward half sat. constant for Ca++ ATPase (mM)
Krb=1.8; % backward half sat. constant for Ca++ ATPase (mM)
KSR=1.0; % scaling factor for Ca++ ATPase
Nfb=0.75; % foward cooperativity constant for Ca++ ATPase
Nrb=0.75; % backward cooperativity constant for Ca++ ATPase
vmaxf=1.53*137.0e-6; % Ca++ ATPase forward rate parameter (mM/ms)
vmaxr=1.53*137.0e-6; % Ca++ ATPase backward rate parameter (mM/ms)
 
fb = power(CaiArray./Kfb,Nfb);
rb = power(CaNSRArray./Krb,Nrb);
Jup = KSR.*(vmaxf.*fb - vmaxr.*rb)./(1.0 + fb + rb); 
%Jup = zeros(1,n_shells);


if(~isreal(Jup))
    pause;
end
% Calculate diffusive fluxes
% alpha_myo = 0.002*diff_factor;
% alpha_NSR = 0.002*diff_factor;
alpha_myo = 0.7;
alpha_NSR = 0.7;
diff_flag = 1;

 if(diff_flag)
  for s = 1: n_shells
    radius = (r(s) + r(s+1)) / 2;
    if(s == 1)
        JMyo_diff(s) = (alpha_myo / (radius* dr* dr))* ( r(s)*( CaiArray(s+1) - CaiArray(s)));
        JNSR_diff(s) = (alpha_NSR / (radius* dr* dr))* ( r(s)*( CaNSRArray(s+1) - CaNSRArray(s)));
    else
        if( s == n_shells)
            JMyo_diff(s) = -(alpha_myo / (radius* dr* dr))* ( r(s -1)*( CaiArray(s) - CaiArray(s -1)));
            JNSR_diff(s) = -(alpha_NSR / (radius* dr* dr))* ( r(s -1)*( CaNSRArray(s) - CaNSRArray(s-1)));
        else
            JMyo_diff(s) = (alpha_myo / (radius* dr* dr))* ( (r(s)*( CaiArray(s+1) - CaiArray(s))) - (r(s -1)*( CaiArray(s) - CaiArray(s -1))));
            JNSR_diff(s) = (alpha_NSR / (radius* dr* dr))* ( (r(s)*( CaNSRArray(s+1) - CaNSRArray(s))) - (r(s -1)*( CaNSRArray(s) - CaNSRArray(s-1))));
        end
    end
  end
 else
     JMyo_diff = zeros(1, n_shells);
     JNSR_diff = zeros(1, n_shells);
 end

a1 = kltrpn_minus.*LTRPNCaArray;
dLTRPNCaArray = kltrpn_plus.*CaiArray.*(ones(1,n_shells) - LTRPNCaArray) - a1;

a1 = khtrpn_minus.*HTRPNCaArray;
dHTRPNCaArray = khtrpn_plus.*CaiArray.*(ones(1,n_shells) - HTRPNCaArray) - a1;
Jtrpn = LTRPNtot.*dLTRPNCaArray + HTRPNtot.*dHTRPNCaArray;

%Jtrpn = zeros(1,n_shells);    
a1 = (CMDNtot*KmCMDN.*ones(1,n_shells))./((CaiArray + KmCMDN.*ones(1,n_shells)).^2);
a2 = (EGTAtot*KmEGTA.*ones(1,n_shells))./((CaiArray +KmEGTA.*ones(1,n_shells)).^2);
beta_i = 1.0./(1.0+a1+a2);

% if(abs(sum(JNSR_diff)) > tol || abs(sum(JMyo_diff)) > tol)
%     disp('Error in diffusion');
% end
dCaiArray = zeros(1,n_shells);
dCaNSRArray = zeros(1, n_shells);
dCaiArray(1) = JMyo_diff(1);
dCaNSRArray(1) = JNSR_diff(1);
for s = 2:n_shells
        dCaiArray(s) = beta_i(s)*(-Jup(s) -Jtrpn(s) + JRyR_int(s-1)*(VCSR/Vmyo(s))) + JMyo_diff(s);
        dCaNSRArray(s) = Jup(s)*(Vmyo(s)/VNSR(s))  - Jtr_int(s-1)*(VCSR/VNSR(s)) + JNSR_diff(s);
        %dCaiArray(s) = JMyo_diff(s);
        %dCaNSRArray(s) = JNSR_diff(s);
end

end

