function [u,ud1,varargout]=eshelbymechfield(a,x,lam,es,nu)

x=col(x);
a=col(a);

%% Determine lambda derivatives

ld1=zeros(3,1);
ld2=zeros(3,3);
Ci=zeros(3,1);
Fij=zeros(3,3);

pa = prod(a);
la=zeros(1,3);
for i = 1 : 3
    la(i) = a(i)^2 + lam;
end
pla = prod(la);
spla=sqrt(pla);

if lam>0
    
    C = sum(x.*x./(a.^2+lam).^2);
    
    for i = 1 : 3
        ld1(i) = 2*x(i)/(la(i))/C;
    end
    
    for i = 1 : 3
%         Ci(i) = x(i)^2/(la(i))^2-sum(x.*x./(a.^2+lam).^3)*2*x(i)/(la(i))/C;
        Ci(i) = 2*x(i)/(la(i))^2-2*sum(x.^2./(a.^2+lam).^3)*2*x(i)/(la(i))/C;
        for j = 1 : 3
            Fij(i,j) = 2*krondelt(i,j)/(la(i))-2*x(i)/(la(i))^2*ld1(j);
        end
    end
    
    
    for i = 1 : 3
        for j = 1 : 3
            ld2(i,j) = (Fij(i,j) - ld1(i)*Ci(j))/C;
        end
    end

end


%% Determine I-integrals

iis=iint(a,lam);
i0=iis{1};
i1=[iis{2}];
i2=[iis{3}];

%% Determine I-integral derivatives

i1d1=zeros(3,3);
i1d2=zeros(3,3,3);
i2d1=zeros(3,3,3);
i2d2=zeros(3,3,3,3);

for i = 1 : 3
    for j = 1 : 3
        i1d1(i,j) = -2*pi*pa/la(i)/spla*ld1(j);
        for k = 1 : 3
            i1d2(i,j,k) = -2*pi*pa/la(i)/spla*(ld2(j,k)-ld1(j)*ld1(k)*(1/la(i)+0.5*sum(1./la)));
            i2d1(i,j,k) = -2*pi*pa/la(i)/la(j)/spla*ld1(k);
            for l = 1 : 3
                i2d2(i,j,k,l) = -2*pi*pa/la(i)/la(j)/spla*(ld2(k,l)-ld1(k)*ld1(l)*(1/la(i)+1/la(j)+0.5*sum(1./la)));
            end
        end
    end
end


%% Determine mechanical fields

psd3=zeros(3,3,3);
psd4=zeros(3,3,3,3);
phd1=zeros(3,1);
phd2=zeros(3,3);

kd=eye(3);

for i = 1 : 3
    phd1(i) = -x(i)*i1(i);
    for j = 1 : 3
        phd2(i,j) = -kd(i,j)*i1(i) - x(i)*i1d1(i,j);
        for l = 1 : 3
            psd3(i,j,l) = -kd(i,j)*x(l)*(i1(l)-a(i)^2*i2(i,l)) - x(i)*x(j)*(i1d1(j,l)-a(i)^2*i2d1(i,j,l)) - (kd(i,l)*x(j)+kd(j,l)*x(i))*(i1(j)-a(i)^2*i2(i,j));
            for k = 1 : 3
                psd4(i,j,k,l) = -kd(i,j)*kd(k,l)*(i1(k)-a(i)^2*i2(i,k))     - (kd(i,k)*kd(j,l)+kd(j,k)*kd(i,l))*(i1(j)-a(i)^2*i2(i,j)) ...
                                -kd(i,j)*x(k)*(i1d1(k,l)-a(i)^2*i2d1(i,k,l)) - (kd(i,k)*x(j)+kd(j,k)*x(i))*(i1d1(j,l)-a(i)^2*i2d1(i,j,l)) ...
                               -(kd(i,l)*x(j)+kd(j,l)*x(i))*(i1d1(j,k)-a(i)^2*i2d1(i,j,k)) - x(i)*x(j)*(i1d2(j,k,l)-a(i)^2*i2d2(i,j,k,l));
            end
        end
    end
end

u=zeros(3,1);
ud1=zeros(3,3);

sdes=sum(diag(es));

for i = 1 : 3
    u(i) = u(i) + 1/8/pi/(1-nu)*(-2*nu*sdes*phd1(i));
    for k = 1 : 3
        u(i) = u(i) + 1/8/pi/(1-nu)*(-4*(1-nu)*es(i,k)*phd1(k));
        ud1(i,k) = 1/8/pi/(1-nu)*(-2*nu*sdes*phd2(i,k));
        for j = 1 : 3
            u(i) = u(i) + 1/8/pi/(1-nu)*(es(k,j)*psd3(k,j,i));
            ud1(i,k) = ud1(i,k) + 1/8/pi/(1-nu)*(-4*(1-nu)*es(i,j)*phd2(j,k));
            for l = 1 : 3
                ud1(i,k) = ud1(i,k) + 1/8/pi/(1-nu)*(es(j,l)*psd4(j,l,i,k));
            end
        end
    end
end

if (nargout > 2)
    varargout{1}=phd1;
    if (nargout > 3)
        varargout{2}=psd3;
        if (nargout > 4)
            varargout{3}=phd2;
            if (nargout > 5)
                varargout{4}=psd4;
            end
        end
    end
end
