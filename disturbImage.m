function Bnew = disturbImage(B,q)
% disturbs image B by flipping pixels with probability q
    [nr nc] = size(B);
    U = unifrnd(0,1,nr,nc);
    Bnew = (U<q).*(1-B)+(U>=q).*B;
end