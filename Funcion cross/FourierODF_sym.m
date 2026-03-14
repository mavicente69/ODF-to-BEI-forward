function [CL,ind_CL]= FourierODF_sym(odf,CS,SS,Lmax)
%% Fourier coeff symetrizados
% odf odf
% Lmax order of the maximum value of L
% CS,SS crystal and sample symmetry
% CL Fourier coeff
% ind_CL order of Foueir coef
% se utiliza la convencion de Bunge para los Tlmn=(-1)^(n+m)Dlmn((alfa,beta,gamma)^-1)

% SO3=get(odf,'center');
% quat=get(SO3,'grid');
% fracc=getdata(odf);
% CS_odf=get(ODF,'CS');
% SS_odf=get(ODF,'CS');
% ker=getkernel(odf);
% 
% quat=quat';  % invierto los quaterniones para expresar los angulos en phi1,Phi,phi2
% odf_inv=ODF(quat,fracc,ker,CS_odf,SS_odf);
% coef_Four=Fourier(odf_inv,'bandwidth',Lmax);

fodf = FourierODF(odf,Lmax);
 
coef_Four=fodf.components{1}.f_hat;

% ind_M=[];
% for l=0:Lmax
%     for m=-l:l
%         for n=-l:l
% ind_M=[ind_M; l m n];
%         end;
%     end;
% end;
% 
% Clmn=[ind_M real(coef_Four) imag(coef_Four)];

% Identificacion del indice donde comienzan los coef con orden l
inicio=1; ind_l=[1];
for i=0:Lmax
   inicio=inicio+(2*i+1)^2;
   ind_l=[ind_l inicio];
end;
% 
% odf_Fourier=FourierODF(coef_Four,CS,SS);
% figure(1001); plotpdf(odf_Fourier,[Miller(1,1,1) Miller(2,0,0) Miller(2,2,0)],'antipodal')
CL=[]; ind_CL=[];
for L=0:Lmax
Ammu=Ammu_sym(L, CS);
Annu=Ammu_sym(L, SS);
if size(Ammu,2)*size(Annu,2)>0
   
c_F=reshape(coef_Four(ind_l(L+1):ind_l(L+1)+(2*L+1)^2-1),(2*L+1),(2*L+1))';

% Transformo entre los Dlmn y los Tlmn mediante la relacion
% Tlmn=(-1)^(n+m)*Dlmn(g^-1). Uso la relacion Tlmn(g^-1)=conj(Tlnm(g))
mm=(sign([-L:L])).^([-L:L])'*(sign([-L:L])).^([-L:L]);
c_F=mm.*c_F;    % multiplico por (-1)^(n+m)
c_F=c_F'; % coef de Fourier !!!! en MTEX5 no se conjuga

CLmunu=Ammu'*c_F*Annu;
CL=[CL; reshape(CLmunu,size(CLmunu,1)*size(CLmunu,2),1)];
for nu=1:size(CLmunu,2)
    for mu=1:size(CLmunu,1)
ind_CL=[ind_CL; L mu nu];
    end; 
end;
    end;
end

end



