function pf=pf_Fourier_sym(CL,ind_CL,CS,SS,iMiller,grid)
% Varios indices de Miller
% Definicion de la direccion hkl

aBR=CS.axes;
v0=dot(aBR(1),cross(aBR(2),aBR(3)));
b1=2*pi/v0*cross(aBR(2),aBR(3));
b2=2*pi/v0*cross(aBR(1),aBR(3));
b3=2*pi/v0*cross(aBR(1),aBR(2));
    
h=iMiller.h; k=iMiller.k; l=iMiller.l;

if strcmp(CS.pointGroup,'6/mmm')
G=h*b1+(h+k)*b2+l*b3;     
else
G=h*b1+k*b2+l*b3; 
end
[rh,th]=polar(G);

[rho,theta] = polar(grid);   % polar, azimuthal
t=reshape(theta,numel(theta),1); r=reshape(rho,numel(rho),1);

PF_hkl=0*ones(length(h),length(t));

for L=0: max(ind_CL(:,1))
j_L=find(ind_CL(:,1)==L);
if numel(j_L)>0
sph_L_mu=sphericalY_sym(L,rh,th,CS);  % radial, azimuthal
sph_L_nu=sphericalY_sym(L,r,t,SS);  % radial, azimuthal
n_mu=size(sph_L_mu,2); n_nu=size(sph_L_nu,2);
CL_M=reshape(CL(j_L),n_mu,n_nu);

PF_hkl=PF_hkl+4*pi/(2*L+1)*conj(sph_L_mu)*CL_M*sph_L_nu';

end
end

CS_pf=crystalSymmetry('triclinic'); CS_pf=CS;
SS_pf=specimenSymmetry('triclinic'); SS_pf=SS;

for jj=1:length(h)
pf{jj}=PoleFigure(iMiller(jj),grid,real(PF_hkl(jj,:)),CS_pf,SS_pf); 
end


pf = [pf{:}];

end