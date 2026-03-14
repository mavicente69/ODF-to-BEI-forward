clear all;


a=4.0479; % Al
%a=3.6151; % Cu

CS=crystalSymmetry('cubic',[a,a,a]); 
SS=specimenSymmetry('triclinic');

h=[1 2 2 3 2 4 4 3];
k=[1 0 2 1 2 0 2 3];
l=[1 0 0 1 2 0 0 1];

x_pos=2*a./sqrt(h.^2+k.^2+l.^2);
y_pos=7*ones(1,length(x_pos));
texto={'111' '200' '220' '311' '222' '400' '420' '331' };
%text(x_pos,y_pos,texto,'Rotation',90)


%% Datos experimentales Nuevos - Genero textura
filename='C:\Users\mavic\MiguelAngel\Modelo de Seccion Eficaz Coherente\MATLAB Simetrizado\Programa Modelo Inversion - MTEX5\Datos Cu nuevo\Miguel_Cu_new.mat';
load(filename);

ker=deLaValleePoussinKernel('HALFWIDTH',5*degree);

pesos=ODF(:,4); 
ori = orientation.byEuler(ODF(:,1),ODF(:,2),ODF(:,3),CS,SS,'Bunge'); clear ODF;

% pesos=1; 
% ori = orientation.byEuler(0,0,0,CS,SS,'Bunge'); clear ODF;

% the model odf
odf = unimodalODF(ori,'weights',pesos,'halfwidth',5*degree);

% lets plot some pole figures
figure(2)
h = [Miller(1,1,1,CS),Miller(2,0,0,CS),Miller(2,2,0,CS)];
plotPDF(odf,h,'antipodal','silent')

% Puntos donde se realizo la medicion del Cu
phi=BB(:,1) ; beta=BB(:,2);
v_exp=vector3d.byPolar(phi,beta);

figure(2); hold on;
id_point={};
for ii=1:length(v_exp)
    id_point=[id_point, num2str(ii)];
end
plot(v_exp,'MarkerEdgeColor','k','MarkerFaceColor','w','Marker','s','label',id_point,'antipodal'); hold off;

%% Pole figure
resol_ang=5*degree;
grid=regularS2Grid('resolution',resol_ang);

[rho,theta] = polar(grid);   % polar, azimuthal
beta=reshape(theta,numel(theta),1); phi=reshape(rho,numel(rho),1);
v_exp=vector3d.byPolar(phi,beta);

%% Calculo de seccion eficaz elastica-coherente
 step_d=0.005; lambda_t_min=1; lambda_t_max=6;
 filename='Al_lattice.txt';
 %filename='Al_latticev2.txt';
 %filename='Cu_lattice.txt';
 %filename='Cu_latticev2.txt';
 tau=0.02; % en Amstrong 
 uhkl=0.02;
 inst_resol=[lambda_t_min 0.01 0.01; 2.5 0.01 0.01; lambda_t_max 0.01 0.01];

[lambda_teor,sigma_teor,m_Y_mat,ind_x,CL,ind_CL]=cross_el_coh(odf,v_exp,filename,lambda_t_min,lambda_t_max,step_d,inst_resol);

%% Cargo dato de seccion eficaz no elastica coherente
filename='C:\Users\mavic\MiguelAngel\Modelo de Seccion Eficaz Coherente\Caso AL7150\Datos EnginX\theory_Powder_cell.dat';
[lambda,isotropico,noelastico]=textread(filename,'%f %f %f','headerlines',1);
figure(16); plot(lambda,[isotropico noelastico])
isotropico_int=interp1(lambda,isotropico,lambda_teor);
noelastico_int=interp1(lambda,noelastico,lambda_teor);
figure(77); plot(lambda_teor,[(isotropico_int-noelastico_int)' sigma_teor(:,1)]); text(x_pos,y_pos,texto,'Rotation',90)


%% Calcula las figuras de olos a partir de la altura de los bordes de Bragg
h=1; k=1; l=1; dhkl=a/sqrt((h^2+k^2+l^2)); i_111=find(abs(lambda_teor-2*dhkl)<step_d/2);
h=2; k=2; l=0; dhkl=a/sqrt((h^2+k^2+l^2)); i_220=find(abs(lambda_teor-2*dhkl)<step_d/2);
h=0; k=0; l=2; dhkl=a/sqrt((h^2+k^2+l^2)); i_002=find(abs(lambda_teor-2*dhkl)<step_d/2);

clear pf_se
pf_se{1}=PoleFigure(Miller(1,1,1,CS),grid,sigma_teor(i_111-2,:)-sigma_teor(i_111+2,:),crystalSymmetry('triclinic'),specimenSymmetry('triclinic'));
pf_se{2}=PoleFigure(Miller(1,1,0,CS),grid,sigma_teor(i_220-2,:)-sigma_teor(i_220+2,:),crystalSymmetry('triclinic'),specimenSymmetry('triclinic'));
pf_se{3}=PoleFigure(Miller(0,0,2,CS),grid,sigma_teor(i_002-2,:)-sigma_teor(i_002+2,:),crystalSymmetry('triclinic'),specimenSymmetry('triclinic'));
pf_se = [pf_se{:}];

figure(9); plot(pf_se)

