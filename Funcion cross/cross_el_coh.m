function [lambda_t,sigma_teor,m_Y_mat,ind_x,CL,ind_CL]=cross_el_coh(odf,v_exp,filename,lambda_t_min,lambda_t_max,step_d,inst_resol)

% This function evaluates the total coherent elastic cross section in textured
% materials in the kinematic approximation in transmission geommetry
% Output
% sigma_teor coherent elastic neutron cross section as a matrix
% m_Y_mat matrix for cross section inversion
% ind_x indeces of m_Y_mat 
% Inputs
% odf Orientation Distribution Function (MTEX class)
% v_exp direction of the beam in the sample refernce system (MTEX vector3d class)
% filename defines the material crystallogrphy and scattering legths
% values of lambda of neutron to be used for calcuation
% tau, uhkl  resolution function of the instrument


CS=odf.CS;
SS=odf.SS;

%% Fourier coeff simetrizados
Lmax=10;
[CL,ind_CL]= FourierODF_sym(odf,CS,SS,Lmax);

% Se utiliza si el MTEX reduce el maximo de bandwidth
Lmax=max(ind_CL(:,1));



%% Definicion de matrix B 
lambda_t=[lambda_t_min:step_d:lambda_t_max];
 % Resolucion Instrumental tipo delta
  if exist('m_Y_mat.mat', 'file') == 0
 [m_Y_mat,ind_x]=matrix_B_sym_Lattice(Lmax,filename,lambda_t);
 save 'm_Y_mat.mat' m_Y_mat ind_x;
 else
 load('m_Y_mat.mat');    
 end


 %% Correccion por funcion  de resolucion instrumental
tau=max(inst_resol(:,2));
uhkl=max(inst_resol(:,3));

rango_l=[-4*(tau+uhkl):step_d:4*(tau+uhkl)]; % rango de valores que definen la funcion instrumental
pos_0=find(abs(rango_l)<step_d/2);
f_Kropff=1/2/tau*exp(-rango_l/tau+uhkl^2/2/tau^2).*erfc(-rango_l/sqrt(2)/uhkl+uhkl/tau); % debo evaluarlo en -lambda por como se hace el calculo
f_resol=f_Kropff'/sum(f_Kropff*step_d);
figure(888); hold on; plot(rango_l,f_resol); hold off;

 m_Y_mat_conv=[]; lambda_t_conv=[];
 
 for j=pos_0:size(m_Y_mat,1)-length(rango_l)+pos_0

     % Dependencia con la long de onda de los tau y uhkl
% tau_l=tau+(lambda_t(j)-4)*0.003;
% uhkl_l=uhkl+(lambda_t(j)-4)*0.003;

tau_l=interp1(inst_resol(:,1),inst_resol(:,2),lambda_t(j));
uhkl_l=interp1(inst_resol(:,1),inst_resol(:,3),lambda_t(j));

f_Kropff=1/2/tau_l*exp(rango_l/tau_l+uhkl_l^2/2/tau_l^2).*erfc(rango_l/sqrt(2)/uhkl_l+uhkl_l/tau_l); % debo evaluarlo en -lambda o lmabda, no se 
f_resol=f_Kropff'/sum(f_Kropff*step_d);
% figure(889); hold on; plot(rango_l,f_resol); 

 m_Y_mat_conv=[m_Y_mat_conv m_Y_mat([j-pos_0+1:j-pos_0+length(rango_l)],:)'*f_resol*step_d]; 
lambda_t_conv=[lambda_t_conv lambda_t(j)];
 end
 
m_Y_mat=m_Y_mat_conv'; lambda_t=lambda_t_conv;
 

%% Seccion eficaz teorica 
% Calculo los coef de la Matriz A teroicos a partir de la textura

[Phi,Beta] = polar(v_exp);   % polar, azimuthal

sigma_teor=[];

for i_g=1:length(Phi)
phi=Phi(i_g); beta=Beta(i_g);

ALnu_teor=[];
for L=0:max(ind_CL(:,1))
i_L=find(ind_CL(:,1)==L);    
klmu=sphericalY_sym(L,phi,beta,SS);
if numel(i_L)>0
n_sol=max(ind_CL(i_L,2));
CL_M=reshape(CL(i_L),n_sol,numel(i_L)/n_sol);    
ALnu=conj(CL_M)*conj(klmu)';
ALnu_teor= [ALnu_teor; ALnu];
end;
end
sigma_teor=[sigma_teor real(m_Y_mat*ALnu_teor)];
end


end


