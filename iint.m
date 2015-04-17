function [iints]=iint(a,lam)

[b ix]=sort(a,'descend');

pa = prod(a);

da12 = a(1)^2 - a(2)^2;
da13 = a(1)^2 - a(3)^2;
da23 = a(2)^2 - a(3)^2;

lv=1e-6;

la=zeros(1,3);
for i = 1 : 3
    la(i) = a(i)^2 + lam;
end

if abs(a(1)-a(2))<lv && abs(a(2)-a(3))<lv % sphere
    
    iints{1} = 4*pi*pa/sqrt(a(1)^2+lam);
    iints{2} = 4*pi*pa/3/(a(1)^2+lam)^1.5*ones(1,3);
    iints{3} = 4*pi*pa/5/la(1)^2.5*ones(3,3);
    
elseif abs(a(1)-a(2))<lv % oblate ellipsoid a1=a2>a3
    
    th = asin(sqrt(da13/la(1)));
    iints{1} = 4*pi*pa/sqrt(da13)*th;
    tv(1) = 2*pi*pa/da13^1.5*(th-sin(2*th)/2);
    tv(2) = tv(1);
    tv(3) = 4*pi*pa/da23/sqrt(da13)*(sqrt(da13)/sqrt(la(3))-th);
    iints{2} = tv;
    
    tm=zeros(3,3);
    tm(1,3) = (tv(3) - tv(1))/(a(1)^2-a(3)^2);
    tm(3,1) = tm(1,3);
    tm(2,3) = tm(1,3);
    tm(3,2) = tm(2,3);
    tm(3,3) = 4*pi/3*pa/la(3)/prod(sqrt(la))-1/3*(tm(3,1)+tm(3,2));
    tm(1,1) = pi*pa/la(1)/prod(sqrt(la))-0.25*tm(1,3);
    tm(1,2) = tm(1,1);
    tm(2,1) = tm(1,2);
    tm(2,2) = tm(1,1);
    
    iints{3} = tm;
        
elseif abs(a(2)-a(3))<lv % prolate ellipsoid a1>a2=a3
    
%     fv = (sqrt(la(1))-sqrt(da13))/(sqrt(la(1))+da13);
%     fv = (sqrt(la(1))+sqrt(la(3)))/(sqrt(la(1))-sqrt(la(3)));
    fv=sqrt(da13)/sqrt(la(3))+sqrt(la(1))/sqrt(la(3));
    lfv=log(fv);
%     lfv=abs(log(fv));
    
    iints{1} = 4*pi*pa/sqrt(da13)*lfv;
    
    tv = zeros(1,3);
    tv(1) = 4*pi*pa/da13^1.5*(lfv-sqrt(da13)/sqrt(la(1)));
    tv(2) = 2*pi*pa/da13^1.5*(sqrt(la(1))*sqrt(da13)/la(3)-lfv);
    tv(3) = tv(2);
    iints{2} = tv;
    
    tm=zeros(3,3);
    tm(1,2) = (tv(2) - tv(1))/(a(1)^2-a(2)^2);
    tm(2,1) = tm(1,2);
    tm(1,3) = tm(1,2);
    tm(3,1) = tm(1,3);
    tm(1,1) = 4*pi/3*pa/la(1)/prod(sqrt(la))-1/3*(tm(1,2)+tm(1,3));
    tm(2,2) = pi*pa/la(2)/prod(sqrt(la))-0.25*tm(1,2);
    tm(2,3) = tm(2,2);
    tm(3,2) = tm(2,3);
    tm(3,3) = tm(2,2);
    
    iints{3} = tm;
    
else % ellipsoid a1>a2>a3
    
    th = asin(sqrt(da13/la(1)));
    k2 = da12/da13;
    
    E = ellipticE(th,k2);
    F = ellipticF(th,k2);
    
    iints{1} = 4*pi*pa/sqrt(da13)*F;
    
    tv = zeros(1,3);
    tv(1) = 4*pi*pa/da12/sqrt(da13)*(F-E);
    tv(2) = 4*pi*pa*(sqrt(da13)/da12/da23*E-1/da12/sqrt(da13)*F-sqrt(la(3))/da23/sqrt(la(1))/sqrt(la(2)));
    tv(3) = 4*pi*pa/da23/sqrt(da13)*(sqrt(la(2))*sqrt(da13)/sqrt(la(1))/sqrt(la(3))-E);
    iints{2} = tv;
    
    tm=zeros(3,3);
    for i = 1 : 3
        for j = 1 : 3
            if i ~= j
                tm(i,j) = (tv(j) - tv(i))/(a(i)^2-a(j)^2);
            end
        end
        tm(i,i) = 4*pi/3*pa/la(i)/prod(sqrt(la))-1/3*(tm(i,mod(i,3)+1)+tm(i,mod(i+1,3)+1));
    end
    iints{3} = tm;
    

end

end