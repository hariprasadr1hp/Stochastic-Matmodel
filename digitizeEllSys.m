function B = digitizeEllSys(XYABP,W,spacing)
% digitizeEllSys returns a 2D matrix
    nrow = ceil((W(2,2)-W(2,1))/spacing);
    ncol = ceil((W(1,2)-W(1,1))/spacing);
    xCoord = repmat(((1:ncol)-0.5)*spacing,nrow,1);
    yCoord = repmat((((1:nrow)-0.5)*spacing)',1,ncol);
    [n nc] = size(XYABP);
    M = XYABP;
    M(:,1) = M(:,1)-W(1,1);
    M(:,2) = M(:,2)-W(2,1);
    B = zeros(nrow,ncol);
    for i=1:n
        Mmax = max(M(i,3),M(i,4));
        imin = max(1,floor((M(i,2)-Mmax+0.5*spacing)/spacing));
        jmin = max(1,floor((M(i,1)-Mmax+0.5*spacing)/spacing));
        imax = min(nrow,floor((M(i,2)+Mmax+0.5*spacing)/spacing));
        jmax = min(ncol,floor((M(i,1)+Mmax+0.5*spacing)/spacing));
        if (imax>=1 && jmax>=1 && imin<=nrow && jmin<=ncol)
            B(imin:imax,jmin:jmax) = max(B(imin:imax,jmin:jmax),((xCoord(imin:imax,jmin:jmax)-M(i,1))*cos(M(i,5))+(yCoord(imin:imax,jmin:jmax)-M(i,2))*sin(M(i,5))).^2/M(i,3)^2+((xCoord(imin:imax,jmin:jmax)-M(i,1))*sin(M(i,5))-(yCoord(imin:imax,jmin:jmax)-M(i,2))*cos(M(i,5))).^2/M(i,4)^2<=1);
        end
    end
end