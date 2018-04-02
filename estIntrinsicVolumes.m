function ALX=estIntrinsicVolumes(B,spacing)
% estIntrinsicVolumes returns a 3 vector 
% containing the 3 quermass densities A_A, L_A, chi_A of B
    F1=[[1 4];[2 8]];
    s=size(B);
    n0=prod(s-1);
    G1=B(1:(s(1)-1),1:(s(2)-1))*F1(1,1)+B(2:s(1),1:(s(2)-1))*F1(2,1)+B(1:(s(1)-1),2:s(2))*F1(1,2)+B(2:s(1),2:s(2))*F1(2,2);

    h=zeros(16,1);
    for k=0:15; h(k+1,1)=sum(sum(G1==k)); end

    theta=pi*[0 0.25 0.5 0.75 1 1.25 1.5 1.75];
    r=[1 sqrt(2) 1  sqrt(2) 1 sqrt(2) 1 sqrt(2)];
    kappa0=[1 1 4 4 8 8 2 2];
    kappa1=[4 8 8 2 2 1 1 4];
    k=1:16;
    pA=zeros(1,8);
    for l=1:8; pA(l)=(bitor(k-1,kappa0(l))==(k-1)).*(bitand(k-1,kappa1(l))==0)*h/r(l)/n0; end

    c=zeros(1,8);
    c(1)=(2*pi-theta(8)+theta(2))/4/pi;
    c(8)=(2*pi-theta(7)+theta(1))/4/pi;
    for l=2:7; c(l)=(theta(l+1)-theta(l-1))/4/pi; end

    AA=(bitor(k-1,1)==(k-1))*h/n0;
    LA=pi*pA*c';
    chiA=(h(2)+h(3)+h(5)+h(9)-h(8)-h(12)-h(14)-h(15))/4/n0;

    ALX=[AA, LA/spacing, chiA/spacing^2];
end