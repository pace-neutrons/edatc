function weight=Bethe_Ansatz_CuGeO3(qh,qk,ql,energy,pars)
%
% Function to calculate the cross-section the Bethe Ansatz for 1d chains,
% e.g. KCuF3, or CuGeO3

%==============
%Muller ansatz code from Tobyfit (Fortran):
% qsqr = qx**2 + qy**2 + qz**2
% wl = pi*p(2)*abs(sin(twopi*ql))
% wu = twopi*p(2)*abs(sin(pi*ql))
% if (eps .gt. (wl+small) .and. eps .le. wu) then
%     weight = p(1) * (bose(eps,temp)/eps) * (form_table(qsqr))**2 / sqrt(eps**2-wl**2)
% else
%     weight = 0.0d0
% endif
%==============

%Avoid rounding issues
small=1e-10;

%Define that we work at base temperature (could set this as a parameter if
%needed):
temp=2;

%Form factor calc
modQ=sqrt(((2*pi/4.81).*qh).^2 + ((2*pi/8.47).*qk).^2 + ((2*pi/2.94).*ql).^2);%Hard wired for CuGeO3
% qsqr=modQ.^2;
FF = Cu_FormFactor(modQ);

%Bethe/Muller Ansatz
wl = pi*pars(2)*abs(sin(2*pi*ql));%choose chains along c
wu = 2*pi*pars(2)*abs(sin(pi*ql));

%weight = pars(1) .* (bose(energy,temp)./energy) .* FF.^2 ./ sqrt(energy.^2-wl.^2);
weight = pars(1)  .* FF.^2 ./ sqrt(energy.^2-wl.^2);
weight(energy<wl+small)=0;
weight(energy>wu)=0;

%Apply a broadening. Need to account for input as a vector or scalar. For
%former assume vector is energy...
% if isvector(weight)
%     gconv=gauss([-10:10],[1,0,5]);
%     weight=conv(weight,gconv,'same');
% else
%     [XX,YY]=meshgrid([-10:10],[-10:10]);
%     gconv=gauss2d(XX,YY,[1,0,0,2,0,2]);
%     weight=conv2(weight,gconv,'same');
% end
end

%==========================================================================
%Tobyfit style Bose factor calculation
%==========================================================================
function BoseFac=bose(energy,temp)

if energy>=0
    BoseFac = energy./(1-exp(-11.6.*energy./temp));
else
    BoseFac = exp(-abs(energy./temp)).*abs(energy)./(1-exp(-11.6.*abs(energy)./temp));
end

end

%==========================================================================
%Cu2+ form factor calculation
%==========================================================================
function FF = Cu_FormFactor(modQ)
%
% Calculate the magnetic form factor for O- ions, using the parameters given
% on the ILL website.
%
% Input |Q| and return the form factor.

% [formFactVal, coeff, S] = sw_mff('Cu2+');%requires SpinW software installed and initialised
% 
% FFj0A=coeff(1); FFj0a= coeff(2); FFj0B= coeff(3); FFj0b= coeff(4); FFj0C= coeff(5); FFj0c= coeff(6); FFj0D=coeff(7);
% 
% s=modQ./(4*pi);
% j0=FFj0A.*exp(-FFj0a.*(s.^2)) + FFj0B.*exp(-FFj0b.*(s.^2)) + FFj0C.*exp(-FFj0c.*(s.^2)) + FFj0D;
% FF=j0;

A=0.0263; a=34.9587; B=0.3668; b=15.9435; C=0.6188; c=5.5935; D=-0.0119;
s=modQ./(4*pi);

j0=A.*exp(-a.*(s.^2)) + B.*exp(-b.*(s.^2)) + C.*exp(-c.*(s.^2)) + D;

FF=j0;
end