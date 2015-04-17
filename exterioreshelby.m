function sol=exterioreshelby(varargin)

if nargin==0
    ltyp=2;
    triax=Inf;
    ef=1e-2;
    np=201;
    plane=12;
    a=[.05 .05 .05];
    b=[.5 .5 .5];
else
    if nargin~=7
        fprintf('Expected 7 inputs, received %d\n',nargin);
        sol=struct();
        return
    end
    a=varargin{1};
    b=varargin{2};
    ef=varargin{3};
    ltyp=varargin{4};
    triax=varargin{5};
    np=varargin{6};
    plane=varargin{7};
end

explot=0;

E=1e3;
nu=0.3;

pcp=paramcp;
c11=pcp(12)+pcp(6)*pcp(11);
c12=pcp(14)+pcp(6)*pcp(13);
c44=pcp(16)+pcp(6)*pcp(15);
E=(c11-c12)*(c11+2*c12)/(c11+c12);
nu=c12/(c11+c12);

C=zeros(3,3,3,3);
Cv=zeros(6,6);
pm=E/(1+nu)/(1-2*nu);
Cv(1,1)=pm*(1-nu);
Cv(2,2)=Cv(1,1);
Cv(3,3)=Cv(1,1);
Cv(1,2)=pm*nu;
Cv(1,3)=Cv(1,2);
Cv(2,3)=Cv(1,2);
Cv(2,1)=Cv(1,2);
Cv(3,1)=Cv(1,2);
Cv(3,2)=Cv(1,2);
Cv(4,4)=pm*(1-2*nu)/2;
Cv(5,5)=Cv(4,4);
Cv(6,6)=Cv(4,4);
C=voigt2(Cv);

e0=zeros(3,3);

if ltyp==1
    e0(1,1)=ef; e0(2,2)=-nu*e0(1,1); e0(3,3)=-nu*e0(1,1); % uniaxial stress
    R=eye(3);
elseif ltyp==2
    e0(1,1)=ef; e0(2,2)=e0(1,1); e0(3,3)=e0(1,1);         % hydrostatic stress
    R=eye(3);
elseif ltyp==3
    e0(1,1)=ef;                                           % 1/3 < chi < Inf
    strainratio=((3*triax-1)*Cv(1,1)-(3*triax+2)*Cv(1,2))/((3*triax+2)*Cv(1,1)-(3*triax-4)*Cv(1,2));
    e0(2,2)=strainratio*e0(1,1); e0(3,3)=strainratio*e0(1,1);
    R=eye(3);
elseif ltyp==4
    e0(1,2)=ef; e0(2,1)=e0(1,2);                          % pure shear
    [R eva]=eig(transpose(eye(3)+e0)*(eye(3)+e0));
    R=eye(3);
    e0=transpose(R)*e0*R;
end

% a=row((eye(3)+e0)*col(a0));

e0v=voigt2(e0);

bx=b(1); by=b(2); bz=b(3);

x1p=linspace(-bx,bx,np);
x2p=linspace(-by,by,np);
x3p=linspace(-bz,bz,np);


if (plane == 12)
    x3p=0;
elseif (plane == 13)
    x2p=0;
elseif (plane == 23)
    x1p=0;
end

[xg yg zg]=meshgrid(x1p,x2p,x3p);
xc=num2cell([reshape(xg,[1,numel(xg)]); reshape(yg,[1,numel(yg)]); reshape(zg,[1,numel(zg)])],1);
as=struct('x',xc);

dst=ceil(numel(xg)/100);

%% Determine lambda

extlam=0;
    
tic;
clc
disp('Lambda Calculations')

lg=zeros(size(xg));
if extlam==1
    th=zeros(size(xg));
    k=zeros(size(xg));
end
% n1g=zeros(size(xg));
% n2g=zeros(size(xg));
% n3g=zeros(size(xg));
asi=0;
for i1 = 1 : size(xg,1)
    for i2 = 1 : size(xg,2)
        for i3 = 1 : size(xg,3)
            if ((xg(i1,i2,i3)/a(1))^2+((yg(i1,i2,i3))/a(2))^2+((zg(i1,i2,i3))/a(3))^2<=1)
                lg(i1,i2,i3)=0;
            else
                lg(i1,i2,i3)=lambdasolve(a,[xg(i1,i2,i3) yg(i1,i2,i3) zg(i1,i2,i3)]);
            end
%             asi=(i1-1)*size(xg,2)*size(xg,3)+(i2-1)*size(xg,2)+i3;
            asi=asi+1;
            
            if (mod(asi-1,dst*10) == 0)
                clc
                strlam=sprintf('Lambda Calculations %0.f s (%3d %%)',toc,(asi-1)/dst);
                fprintf('%s\n',strlam);
                fprintf('%.0f s remaining\n',(numel(xg)-asi)*toc/asi);
            end
            
%             n1g(i1,i2,i3)=xg(i1,i2,i3)/(a(1)+lg(i1,i2,i3));
%             n2g(i1,i2,i3)=yg(i1,i2,i3)/(a(2)+lg(i1,i2,i3));
%             n3g(i1,i2,i3)=zg(i1,i2,i3)/(a(3)+lg(i1,i2,i3));
%             ng=[n1g(i1,i2,i3) n2g(i1,i2,i3) n3g(i1,i2,i3)];
%             nng=norm(ng);
%             n1g(i1,i2,i3)=n1g(i1,i2,i3)/nng;
%             n2g(i1,i2,i3)=n2g(i1,i2,i3)/nng;
%             n3g(i1,i2,i3)=n3g(i1,i2,i3)/nng;
            if extlam==1
                th(i1,i2,i3)=asin(sqrt((a(1)^2-a(3)^2)/(a(1)^2+lg(i1,i2,i3))));
                k(i1,i2,i3)=sqrt((a(1)^2-a(2)^2)/(a(1)^2-a(3)^2));
            end
        end
    end
end
if extlam==1
    thc=num2cell(reshape(th,[1,numel(th)]));
    [as.th]=thc{:};
    kc=num2cell(reshape(k,[1,numel(k)]));
    [as.k]=kc{:};
end
lc=num2cell(reshape(lg,[1,numel(lg)]));
[as.l]=lc{:};
% nc=num2cell([reshape(n1g,[1,numel(n1g)]); reshape(n2g,[1,numel(n2g)]); reshape(n3g,[1,numel(n3g)])],1);
% [as.n]=nc{:};

tlam=toc;
strlam=sprintf('Lambda Calculations %.0f s',tlam);
clc;
fprintf('%s\n',strlam);


%% Determine I-integrals

iis=iint(a,0);
i0=iis{1};
i1=[iis{2}];
i2=[iis{3}];

%% Determine eigenstrain

tic;
disp('Eigenstrain Calculations')

es=zeros(3,3);
esv=zeros(6,1);

A=zeros(3,3,3,3);
Av=zeros(6,6);
Ainv=zeros(3,3,3,3);
Avinv=zeros(6,6);
S=zeros(3,3,3,3);
Sv=zeros(6,6);


for i = 1 : 3
    S(i,i,i,i)=3/8/pi/(1-nu)*a(i)^2*i2(i,i)+(1-2*nu)/8/pi/(1-nu)*i1(i);
    if i ~= 1
        S(1,1,i,i)=1/8/pi/(1-nu)*a(i)^2*i2(1,i)-(1-2*nu)/8/pi/(1-nu)*i1(1);
        S(1,i,1,i)=(a(1)^2+a(i)^2)/16/pi/(1-nu)*i2(1,i)+(1-2*nu)/16/pi/(1-nu)*(i1(1)+i1(i));
        S(1,i,i,1)=S(1,i,1,i);
    end
    if i ~= 2
        S(2,2,i,i)=1/8/pi/(1-nu)*a(i)^2*i2(2,i)-(1-2*nu)/8/pi/(1-nu)*i1(2);
        S(2,i,2,i)=(a(2)^2+a(i)^2)/16/pi/(1-nu)*i2(2,i)+(1-2*nu)/16/pi/(1-nu)*(i1(2)+i1(i));
        S(2,i,i,2)=S(2,i,2,i);
    end
    if i ~= 3
        S(3,3,i,i)=1/8/pi/(1-nu)*a(i)^2*i2(3,i)-(1-2*nu)/8/pi/(1-nu)*i1(3);
        S(3,i,3,i)=(a(3)^2+a(i)^2)/16/pi/(1-nu)*i2(3,i)+(1-2*nu)/16/pi/(1-nu)*(i1(3)+i1(i));
        S(3,i,i,3)=S(3,i,3,i);
    end
end
Sv=voigt2(S);

Av=inv(Cv*(eye(6)-Sv))*Cv;
A=voigt2(Av);

esv=Av*e0v;
es=voigt2(esv);

[c0 Lv Mv Sv2]=eshelby2(a,2*[bx;by;bz],voigt2(Cv),16,16);

teig=toc;
streig=sprintf('Eigenstrain Calculations %.0f s',teig);
clc;
fprintf('%s\n',strlam);
fprintf('%s\n',streig);

%% Determine mechanical fields

exout=0;

clc
fprintf('%s\n',strlam);
disp('Mechanical Field Calculations')

tic;

for i = 1 : numel(xg)

    if exout==0
        [u ud1]=eshelbymechfield(a,as(i).x,as(i).l,es,nu);
    else
        [u ud1 phd1 psd3 phd2 psd4]=eshelbymechfield(a,as(i).x,as(i).l,es,nu);
        as(i).phd2=phd2;
        as(i).psd3=psd3;
        as(i).psd4=psd4;
    end
    as(i).u=u;
    as(i).ud1=ud1;
    
    if (mod(i-1,dst) == 0)
        clc
        fprintf('%s\n',strlam);
        fprintf('%s\n',streig);
        strmech=sprintf('Mechanical Field Calculations %0.f s (%3d %%)',toc,(i-1)/dst);
        fprintf('%s\n',strmech);
        fprintf('%.0f s remaining\n',(numel(xg)-i)*toc/i);
    end
    
end
strmech=sprintf('Mechanical Field Calculations %.0f s (%.2g s per dof)',toc,toc/numel(xg));
clc
fprintf('%s\n',strlam);
fprintf('%s\n',streig);
fprintf('%s\n',strmech);

%% Determine moments



%% postprocess

xv=[as.x];
usv=[as.u];
lr=[as.l];
if extlam==1
    kr=[as.k];
end

disp('Postprocessing')
tic

phd2=zeros(6,numel(xg));
epssv=zeros(6,numel(xg));
gamsv=zeros(6,numel(xg));
sigsv=zeros(6,numel(xg));

epshv=repmat(e0v,[1,numel(xg)]);
sig0v=Cv*e0v;
sighv=repmat(sig0v,[1,numel(xg)]).*repmat((sum((xv./repmat(col(a),[1,size(xv,2)])).^2)>1),[6,1]);

for i = 1 : numel(xg)
    if exout~=0
        phd2(:,i)=voigt([as(i).phd2]);
        psd3111(1,i)=as(i).psd3(1,1,1);
        psd41111(1,i)=as(i).psd4(1,1,1,1);
    end
    epssm=0.5*([as(i).ud1]+transpose([as(i).ud1]));
    epssv(:,i)=voigt(epssm);
    uhv(:,i)=e0*xv(:,i);
    if sum((xv(:,i)./col(a)).^2)>1
        gamsv=voigt2(epssm);
        sigsv(:,i)=Cv*gamsv;
    end
end
sigtv=sighv+sigsv; 
epstv=epshv+epssv; 

x1=xv(1,:); x2=xv(2,:); x3=xv(3,:); 

th=atan2(x2,x1); u1=usv(1,:); u2=usv(2,:); u3=usv(3,:); 

if exout~=0
    ph11=phd2(1,:); ph22=phd2(2,:); ph33=phd2(3,:); ph23=phd2(4,:); ph13=phd2(5,:); ph12=phd2(6,:);
end

es11=epssv(1,:); es22=epssv(2,:); es33=epssv(3,:); es23=epssv(4,:);  es13=epssv(5,:); es12=epssv(6,:);  
cs11=sigtv(1,:); cs22=sigtv(2,:); cs33=sigtv(3,:);cs12=sigtv(6,:);cs13=sigtv(5,:);cs23=sigtv(4,:);

cshyd=1/3*(cs11+cs22+cs33);

ds11=cs11-cshyd; ds22=cs22-cshyd; ds33=cs33-cshyd;

vms=sqrt(3/2)*sqrt(ds11.^2+ds22.^2+ds33.^2+2*(cs12.^2+cs13.^2+cs23.^2));

upv=zeros(3,numel(xg));
uet=usv+uhv;

uhr=zeros(1,numel(xg));
usr=zeros(1,numel(xg));
upr=zeros(1,numel(xg));
uer=zeros(1,numel(xg));

epssrr=zeros(1,numel(xg));
sigsrr=zeros(1,numel(xg));
esrr=zeros(1,numel(xg));
fd11=zeros(1,numel(xg));
fd22=zeros(1,numel(xg));
fd33=zeros(1,numel(xg));
fd23=zeros(1,numel(xg));
fd13=zeros(1,numel(xg));
fd12=zeros(1,numel(xg));
sigp11=zeros(1,numel(xg));
sigp22=zeros(1,numel(xg));
sigp33=zeros(1,numel(xg));
sigp23=zeros(1,numel(xg));
sigp13=zeros(1,numel(xg));
sigp12=zeros(1,numel(xg));
sigprr=zeros(1,numel(xg));

for i = 1 : numel(xg)
    
    Rmat=rotmat([0;0;1],-th(i));
    RmatT=transpose(Rmat);
    
    sigsrm=Rmat*voigt(sigsv(:,i))*RmatT;
    sigsrr(i)=sigsrm(1,1);
    
    sighrm=Rmat*voigt(sighv(:,i))*RmatT;
    sighrr(i)=sighrm(1,1);
    
    epssr=Rmat*voigt(epssv(:,i))*RmatT;
    epssrr(i)=epssr(1,1);
    
    epshrm=Rmat*e0*RmatT;
    epshrr(i)=epshrm(1,1);
    
    sdm=[ds11(i) cs12(i) cs13(i); cs12(i) ds22(i) cs23(i); cs13(i) cs23(i) ds33(i)];
    nsd=sqrt(sum(sum(sdm.*sdm)));
    if (nsd > 0)
        fd=sdm/sqrt(sum(sum(sdm.*sdm)));
    else
        fd=zeros(3);
    end
    
    fd11(i)=fd(1,1); fd22(i)=fd(2,2); fd33(i)=fd(3,3); fd23(i)=fd(2,3); fd13(i)=fd(1,3); fd12(i)=fd(1,2);
    
    sigp11(i)=vms(i)*fd(1,1); sigp22(i)=vms(i)*fd(2,2); sigp33(i)=vms(i)*fd(3,3); sigp23(i)=vms(i)*fd(2,3); sigp13(i)=vms(i)*fd(1,3); sigp12(i)=vms(i)*fd(1,2);
    
    sigpm=vms(i)*fd;
    sigpmr=Rmat*sigpm*RmatT;
    sigprr(i)=sigpmr(1,1);
    sigppmr=[sigprr(i) 0 0; 0 0 0; 0 0 0];
    sigppm=RmatT*sigppmr*Rmat;
    
    
    gamp=inv(Cv)*[sigp11(i);sigp22(i);sigp33(i);2*sigp23(i);2*sigp13(i);2*sigp12(i)];
    
    epsp11(i)=gamp(1); epsp22(i)=gamp(2); epsp33(i)=gamp(3); epsp23(i)=0.5*gamp(4); epsp13(i)=0.5*gamp(5); epsp12(i)=0.5*gamp(6);

    epm=voigt2(gamp);
    upv(:,i)=(epm*xv(:,i));
    
    if (norm(xv(:,i)) > 0)
        nr=norm(xv(:,i));
        xvn=xv(:,i)/nr;
    else
        xvn=xv(:,i);
    end
    uhr(i)=sum(uhv(:,i).*xvn);
    usr(i)=sum(usv(:,i).*xvn);
    upr(i)=sum(upv(:,i).*xvn);
    uer(i)=sum(uet(:,i).*xvn);
end

sigpv=[sigp11;sigp22;sigp33;sigp23;sigp13;sigp12];

epstrr = epssrr + epshrr;
sigtrr = sigsrr + sighrr;

nbin=500; 
dl=max(lr)-min(lr); 
lrn=lr/dl; 
lrn=lr/dl+1/(nbin-1)-1e-8*dl;
lrb=floor(lrn*(nbin-1))+1;
ep11b=zeros(size(lr));
ep22b=zeros(size(lr));
ep33b=zeros(size(lr));
ep23b=zeros(size(lr));
ep13b=zeros(size(lr));
ep12b=zeros(size(lr));
for i = 1 : nbin
    binvec=find(lrb==i);
    numbin=numel(binvec);
    ep11av=sum(epsp11(binvec))/numbin; ep11avvec=zeros(size(lr)); ep11avvec(binvec)=ep11av; ep11b=ep11b+ep11avvec;
    ep22av=sum(epsp22(binvec))/numbin; ep22avvec=zeros(size(lr)); ep22avvec(binvec)=ep22av; ep22b=ep22b+ep22avvec;
    ep33av=sum(epsp33(binvec))/numbin; ep33avvec=zeros(size(lr)); ep33avvec(binvec)=ep33av; ep33b=ep33b+ep33avvec;
    ep23av=sum(epsp23(binvec))/numbin; ep23avvec=zeros(size(lr)); ep23avvec(binvec)=ep23av; ep23b=ep23b+ep23avvec;
    ep13av=sum(epsp13(binvec))/numbin; ep13avvec=zeros(size(lr)); ep13avvec(binvec)=ep13av; ep13b=ep13b+ep13avvec;
    ep12av=sum(epsp12(binvec))/numbin; ep12avvec=zeros(size(lr)); ep12avvec(binvec)=ep12av; ep12b=ep12b+ep12avvec;
end

if (plane == 12)
    ivs=find(x3==0);
elseif (plane == 13)
    ivs=find(x2==0);
elseif (plane == 23)
    ivs=find(x1==0);
else
    ivs=find(x3==0);
end

sol.np=np;
sol.plane=plane;
sol.triax=triax;

sol.e0=e0;
sol.xv=xv;
sol.epsp=[epsp11;epsp22;epsp33;epsp23;epsp13;epsp12];

sol.sighv=sighv;
sol.sigsv=sigsv;
sol.sigpv=sigpv;
sol.epshv=epshv;
sol.epssv=epssv;
sol.uhv=uhv;
sol.usv=usv;
sol.upv=upv;
sol.uer=uer;
sol.upr=upr;
sol.sighrr=sighrr;
sol.sigsrr=sigsrr;
sol.sigprr=sigprr;
sol.epssrr=epssrr;
sol.epshrr=epshrr;

strpost=sprintf('Postprocessing %.0f s',toc);
clc
fprintf('%s\n',strlam);
fprintf('%s\n',streig);
fprintf('%s\n',strmech);
fprintf('%s\n',strpost);

%% plot

disp('Plotting')
tic

figure(1)
clf
imagesc(reshape(vms(ivs)/sum(Cv(1:3,:)*e0v)*3,[np,np])); axis equal; axis off
clf
imagesc(reshape(vms(ivs),[np,np])); axis equal; axis off; caxis(max(abs(vms))*[-1 1])
set(gcf,'colormap',[[[0:31]';31*ones(32,1)] [[0:31]';[31:-1:0]'] [31*ones(32,1);[31:-1:0]']]/31)
colorbar;
title(['von Mises Stress (\chi=' num2str(triax,'%.2f') ')'])
set(gca,'fontsize',18)

if explot==1

figure(2)
clf
imagesc(reshape(sqrt(u1(ivs).^2+u2(ivs).^2),[np,np])); axis equal; axis off
colorbar;
title('Magnitude of In-Plane Displacement')

sp=1;
if sp==1
    figure(4); close(4);
    figure(5); close(5);
    figure(6); close(6);
    figure(7); close(7);
    figure(8); close(8);
    figure(10); close(10);
    figure(11); close(11);
    figure(12); close(12);
    figure(13); close(13);
    figure(14); close(14);
end

figure(3)
clf
if sp==1
    subplot(2,3,1)
end
imagesc([-bx bx],[by -by],e0(1,1)+reshape(es11(ivs),[np,np])); axis equal; axis off; %caxis([0.0032093 0.01419]); 
set(gcf,'colormap',[[zeros(21,1);linspace(0,1,22)';ones(21,1)],[linspace(0,1,22)';ones(20,1);linspace(1,0,22)'],[ones(21,1);linspace(1,0,22)';zeros(21,1)]])
colorbar
set(gca,'fontsize',24)
title('Strain \epsilon_{11}')

if sp==1
    subplot(2,3,2)
else
    figure(4)
    clf
end
imagesc([-bx bx],[by -by],e0(2,2)+reshape(es22(ivs),[np,np])); axis equal; axis off; %caxis([-0.00429717 -0.00178362]); 
set(gcf,'colormap',[[zeros(21,1);linspace(0,1,22)';ones(21,1)],[linspace(0,1,22)';ones(20,1);linspace(1,0,22)'],[ones(21,1);linspace(1,0,22)';zeros(21,1)]])
colorbar
set(gca,'fontsize',24)
title('Strain \epsilon_{22}')

if sp==1
    subplot(2,3,3)
else
    figure(5)
    clf
end
imagesc([-bx bx],[by -by],e0(3,3)+reshape(es33(ivs),[np,np])); axis equal; axis off; %caxis([-0.00429717 -0.00178362]); 
set(gcf,'colormap',[[zeros(21,1);linspace(0,1,22)';ones(21,1)],[linspace(0,1,22)';ones(20,1);linspace(1,0,22)'],[ones(21,1);linspace(1,0,22)';zeros(21,1)]])
colorbar
set(gca,'fontsize',24)
title('Strain \epsilon_{33}')

if sp==1
    subplot(2,3,4)
else
    figure(6)
    clf
end
imagesc([-bx bx],[by -by],e0(2,3)+reshape(es23(ivs),[np,np])); axis equal; axis off; %caxis([-0.00603 0.00603]); 
set(gcf,'colormap',[[zeros(21,1);linspace(0,1,22)';ones(21,1)],[linspace(0,1,22)';ones(20,1);linspace(1,0,22)'],[ones(21,1);linspace(1,0,22)';zeros(21,1)]])
colorbar
set(gca,'fontsize',24)
title('Strain \epsilon_{23}')

if sp==1
    subplot(2,3,5)
else
    figure(7)
    clf
end
imagesc([-bx bx],[by -by],e0(1,3)+reshape(es13(ivs),[np,np])); axis equal; axis off; %caxis([-0.00603 0.00603]); 
set(gcf,'colormap',[[zeros(21,1);linspace(0,1,22)';ones(21,1)],[linspace(0,1,22)';ones(20,1);linspace(1,0,22)'],[ones(21,1);linspace(1,0,22)';zeros(21,1)]])
colorbar
set(gca,'fontsize',24)
title('Strain \epsilon_{13}')

if sp==1
    subplot(2,3,6)
else
    figure(8)
    clf
end
imagesc([-bx bx],[by -by],e0(1,2)+reshape(es12(ivs),[np,np])); axis equal; axis off; %caxis([-0.00603 0.00603]); 
set(gcf,'colormap',[[zeros(21,1);linspace(0,1,22)';ones(21,1)],[linspace(0,1,22)';ones(20,1);linspace(1,0,22)'],[ones(21,1);linspace(1,0,22)';zeros(21,1)]])
colorbar
set(gca,'fontsize',24)
title('Strain \epsilon_{12}')

figure(9)
clf
if sp==1
    subplot(2,3,1)
end
imagesc([-bx bx],[by -by],reshape(cs11(ivs),[np,np])); axis equal; axis off;% caxis([0.0032093 0.01419]); 
set(gcf,'colormap',[[zeros(21,1);linspace(0,1,22)';ones(21,1)],[linspace(0,1,22)';ones(20,1);linspace(1,0,22)'],[ones(21,1);linspace(1,0,22)';zeros(21,1)]])
colorbar
set(gca,'fontsize',24)
title('Stress \sigma_{11}')

if sp==1
    subplot(2,3,2)
else
    figure(10)
    clf
end
imagesc([-bx bx],[by -by],reshape(cs22(ivs),[np,np])); axis equal; axis off;% caxis([0 1.325]); 
set(gcf,'colormap',[[zeros(21,1);linspace(0,1,22)';ones(21,1)],[linspace(0,1,22)';ones(20,1);linspace(1,0,22)'],[ones(21,1);linspace(1,0,22)';zeros(21,1)]])
colorbar
set(gca,'fontsize',24)
title('Stress \sigma_{22}')

if sp==1
    subplot(2,3,3)
else
    figure(11)
    clf
end
imagesc([-bx bx],[by -by],reshape(cs33(ivs),[np,np])); axis equal; axis off;% caxis([0 1.325]);
set(gcf,'colormap',[[zeros(21,1);linspace(0,1,22)';ones(21,1)],[linspace(0,1,22)';ones(20,1);linspace(1,0,22)'],[ones(21,1);linspace(1,0,22)';zeros(21,1)]])
colorbar
set(gca,'fontsize',24)
title('Stress \sigma_{33}')

if sp==1
    subplot(2,3,4)
else
    figure(12)
    clf
end
imagesc([-bx bx],[by -by],reshape(cs23(ivs),[np,np])); axis equal; axis off;
set(gcf,'colormap',[[zeros(21,1);linspace(0,1,22)';ones(21,1)],[linspace(0,1,22)';ones(20,1);linspace(1,0,22)'],[ones(21,1);linspace(1,0,22)';zeros(21,1)]])
colorbar
set(gca,'fontsize',24)
title('Stress \sigma_{23}')

if sp==1
    subplot(2,3,5)
else
    figure(13)
    clf
end
imagesc([-bx bx],[by -by],reshape(cs13(ivs),[np,np])); axis equal; axis off;% caxis([-2.346 2.346]); 
set(gcf,'colormap',[[zeros(21,1);linspace(0,1,22)';ones(21,1)],[linspace(0,1,22)';ones(20,1);linspace(1,0,22)'],[ones(21,1);linspace(1,0,22)';zeros(21,1)]])
colorbar
set(gca,'fontsize',24)
title('Stress \sigma_{13}')

if sp==1
    subplot(2,3,6)
else
    figure(14)
    clf
end
imagesc([-bx bx],[by -by],reshape(cs12(ivs),[np,np])); axis equal; axis off;% caxis([-2.346 2.346]); 
set(gcf,'colormap',[[zeros(21,1);linspace(0,1,22)';ones(21,1)],[linspace(0,1,22)';ones(20,1);linspace(1,0,22)'],[ones(21,1);linspace(1,0,22)';zeros(21,1)]])
colorbar
set(gca,'fontsize',24)
title('Stress \sigma_{12}')



figure(15)
clf

if sp==1
    subplot(1,3,1)
end
imagesc([-bx bx],[by -by],reshape(u1(ivs),[np,np])); axis equal; axis off;% caxis([0.0032093 0.01419]); 
set(gcf,'colormap',[[zeros(21,1);linspace(0,1,22)';ones(21,1)],[linspace(0,1,22)';ones(20,1);linspace(1,0,22)'],[ones(21,1);linspace(1,0,22)';zeros(21,1)]])
colorbar
set(gca,'fontsize',24)
title('Displacement u_1')

if sp==1
    subplot(1,3,2)
else
    figure(16)
    clf
end
imagesc([-bx bx],[by -by],reshape(u2(ivs),[np,np])); axis equal; axis off;% caxis([0 1.325]); 
set(gcf,'colormap',[[zeros(21,1);linspace(0,1,22)';ones(21,1)],[linspace(0,1,22)';ones(20,1);linspace(1,0,22)'],[ones(21,1);linspace(1,0,22)';zeros(21,1)]])
colorbar
set(gca,'fontsize',24)
title('Displacement u_2')

if sp==1
    subplot(1,3,3)
else
    figure(17)
    clf
end
imagesc([-bx bx],[by -by],reshape(u3(ivs),[np,np])); axis equal; axis off;% caxis([0 1.325]);
set(gcf,'colormap',[[zeros(21,1);linspace(0,1,22)';ones(21,1)],[linspace(0,1,22)';ones(20,1);linspace(1,0,22)'],[ones(21,1);linspace(1,0,22)';zeros(21,1)]])
colorbar
set(gca,'fontsize',24)
title('Displacement u_3')

figure(20); 
% clf; 
set(gcf,'colormap',[[[0:31]';31*ones(32,1)] [[0:31]';[31:-1:0]'] [31*ones(32,1);[31:-1:0]']]/31)
imagesc([-bx bx],[by -by],reshape(uer,[np,np])); 
axis equal;
caxis(max(max(abs(uer)),max(abs(upr)))*[-1 1])
colorbar; 
set(gca,'fontsize',18)
title('Radial Elastic Displacement u^e_r')

figure(21); 
% clf; 
set(gcf,'colormap',[[[0:31]';31*ones(32,1)] [[0:31]';[31:-1:0]'] [31*ones(32,1);[31:-1:0]']]/31)
imagesc([-bx bx],[by -by],reshape(upr,[np,np]));  
axis equal;
caxis(max(max(abs(uer)),max(abs(upr)))*[-1 1])
colorbar; 
set(gca,'fontsize',18)
title('Radial Plastic Displacement u^p_r')

figure(22); 
% clf; 
set(gcf,'colormap',[[[0:31]';31*ones(32,1)] [[0:31]';[31:-1:0]'] [31*ones(32,1);[31:-1:0]']]/31)
imagesc([-bx bx],[by -by],reshape(uer+upr,[np,np]));  
axis equal;
caxis(max(max(abs(uer)),max(abs(upr)))*[-1 1])
colorbar; 
set(gca,'fontsize',18)
title('Radial Total Displacement u_r')

figure(26); clf; imagesc([-bx bx],[by -by],reshape(epstrr,[np,np])); 
axis equal;
colorbar; 
set(gca,'fontsize',18)
title('Radial Plastic Strain \epsilon^p_{rr}')

figure(27); clf; imagesc([-bx bx],[by -by],reshape(sighrr+sigsrr,[np,np])); 
axis equal;
colorbar; 
set(gca,'fontsize',18)
title('Radial Stress \sigma_{rr}')

figure(28); clf; imagesc([-bx bx],[by -by],reshape(sighrr+sigsrr-sigprr,[np,np])); axis equal; colorbar
set(gca,'fontsize',18)
title('Radial Stress \sigma^e_{rr}-\sigma^p_{rr}')

end

strplot=sprintf('Plotting %.0f s',toc);
clc
fprintf('%s\n',strlam);
fprintf('%s\n',streig);
fprintf('%s\n',strmech);
fprintf('%s\n',strpost);
fprintf('%s\n',strplot);

