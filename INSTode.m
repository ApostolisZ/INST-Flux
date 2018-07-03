function dXdt = INSTode(t, X, INSTmodel)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Z = zeros(size(X));
[Bi, Bj, Bs] = find(INSTmodel.B);
for i = 1:length(Bi)   
    total_mid = 0;
    for isol = 1:size(INSTmodel.Y{Bj(i)}{1}, 1)
        mid = [];
        for jsol = 1:size(INSTmodel.Y{Bj(i)}{1}, 2)
            sxint = deval(INSTmodel.Y{Bj(i)}{1}(isol, jsol), t, INSTmodel.Y{Bj(i)}{2}{isol, jsol});
            mid = conv([ 1 - sum(mid); mid ], [ 1 - sum(sxint); sxint ]);
            mid = mid(2:end);
        end
        total_mid = total_mid + mid;
    end
    
    Z((Bi(i) - 1) + (1:length(total_mid))) = Z((Bi(i) - 1) + (1:length(total_mid))) ...
        + Bs(i) * total_mid;
end

dXdt = (INSTmodel.A * X + Z) ./ INSTmodel.C_diag;
end