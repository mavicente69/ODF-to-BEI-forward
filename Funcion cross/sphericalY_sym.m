function [sph]=sphericalY_sym(L,phi,beta,SS)
% phi angulo polar
% beta angulo azimuthal
sph=[];

[Ammu]=Ammu_sym(L,SS);
if size(Ammu,2)>0
mm=(ones(length(phi),1)*(sign([-L:L]).^([-L:L])));
ylm=sphericalY(L,phi,beta).*mm; % le tuve que restar pi, chequear
sph=ylm*Ammu; % seleccono solo los indep    
end;


end

