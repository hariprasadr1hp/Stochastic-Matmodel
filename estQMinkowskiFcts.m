function ALXq = estQMinkowskiFcts(B,k,spacing)
% estQMinkowskiFcts returns a (k+1) by 4 matrix
% containing the three Minkowski functions in 2D
% according to a square structuring element.
% The first row contains the arguments r.
    s = size(B);
    W = B;
    if k>min(s-3)
        k = min(s-3);
    end
    ALXq = [spacing*(0:k).' zeros(k+1,3)];
    ALXq(1,2:4) = estIntrinsicVolumes(W,spacing);
    if k>0       
        for r=2:(k+1)
            BS = W(1:(s(1)-1),1:(s(2)-1))+W(2:s(1),1:(s(2)-1))+W(1:(s(1)-1),2:s(2))+W(2:s(1),2:s(2));
            W = BS>0;
            ALXq(r,2:4) = estIntrinsicVolumes(W,spacing);
            s = s-1;
        end
    end
end




