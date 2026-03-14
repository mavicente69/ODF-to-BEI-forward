function [m_Y_mat,ind_x]=matrix_B_sym_Lattice(Lmax,filename,lambda_t)

fileID=fopen(filename);
% Lee el archivo que define la red y sus bases
InputText=textscan(fileID,'%s',1,'delimiter','\n'); % Read header line
CS=crystalSymmetry(InputText{1,1});
InputText=textscan(fileID,'%s',1,'delimiter','\n'); % Read header line
RD=textscan(fileID,'%f %f %f',3,'delimiter','\n'); % Read header line
InputText=textscan(fileID,'%s',1,'delimiter','\n'); % Read header line
BS=textscan(fileID,'%f %f %f %f %f %f','delimiter','\n'); % Read header line
fclose(fileID);

% Red real
M=cell2mat(RD); v0=abs(det(M));
a1=vector3d(M(1,1),M(1,2),M(1,3));
a2=vector3d(M(2,1),M(2,2),M(2,3));
a3=vector3d(M(3,1),M(3,2),M(3,3));

% Red Reciproca
b1=2*pi/v0*cross(a2,a3);
b2=2*pi/v0*cross(a3,a1);
b3=2*pi/v0*cross(a1,a2);

% Bases
M=cell2mat(BS(2:4));
r_b=vector3d(M(:,1),M(:,2),M(:,3));
b_scatt=cell2mat(BS(5)); u02=cell2mat(BS(6));

%% Matriz para invertir coeficientes Clmn a ppartir de la seccion eficaz neutronica
 N=1;  % numero de celdas unidades

m_Y_mat=[];


% vectores de la red reciproca sobre lois cuales se hara el calculo
theta_G=[]; rho_G=[]; mod_G=[]; Fhkl2_t=[];
for h=[-40:40];
   for k=[-40:40];
        for l=[-40:40];
             
G=h*b1+k*b2+l*b3; [theta,rho] = polar(G); 
 if norm(G)~=0 && norm(G)<20*norm(b3)
    theta_G=[theta_G; theta]; rho_G=[rho_G; rho]; mod_G=[mod_G; norm(G)];

    fi2=abs(b_scatt'*(exp(complex(0,dot(G,r_b))).*exp(-norm(G)^2*u02/2)))^2/100;
    Fhkl2_t=[Fhkl2_t fi2];
 end
 end
        end
end
sphY_t=[]; ind_x=[];
for L=0:Lmax
sphY=sphericalY_sym(L,theta_G,rho_G,CS);
sphY_t=[sphY_t sphY];
for j=1:size(sphY,2)
ind_x=[ind_x; L j];
end
end
sphY_t=sphY_t';

m_Y_mat=[];
%calculo para distintos lambda   
for lambda=lambda_t
    
ki=2*pi/lambda;
lambda;

% i_G=find(heaviside(1-mod_G/2/ki));
i_G=find(abs(mod_G/2/ki)<=1);

if numel(i_G)>0
MI=sphY_t(:,i_G)*diag(Fhkl2_t(i_G));

Y_t=[];
for i=1:size(ind_x,1)
leg_as=legendre(ind_x(i,1),mod_G(i_G)/2/ki)';
f_mod=(4*pi/(2*ind_x(i,1)+1)*N*(2*pi)^3/v0/ki^3*ki/2).*leg_as(:,1)./mod_G(i_G);
Y_t=[Y_t MI(i,:)*f_mod];
end;
       
else
Y_t=0*ones(1,size(ind_x,1));
end

m_Y_mat=[m_Y_mat; Y_t];


%save 'm_Y_mat_Zr.mat' m_Y_mat ind_x
end

