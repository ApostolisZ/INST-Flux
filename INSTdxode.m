function ddXdt = INSTdxode(t, dX, INSTmodel)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Z = zeros(size(dX));
A = INSTmodel.A(:, cumsum([1; INSTmodel.emu_natoms(1:(end - 1))]));
[Ai, Aj, As] = find(A);
for i = 1:length(Ai)   
    mid = [];
    for isol = 1:length(INSTmodel.X{Aj(i)}{1})
        sxint = deval(INSTmodel.X{Aj(i)}{1}(isol), t, INSTmodel.X{Aj(i)}{2}{isol});
        mid = conv([ 1 - sum(mid); mid ], [ 1 - sum(sxint); sxint ]);
        mid = mid(2:end);
    end
    
    Z((Ai(i) - 1) + (1:length(mid))) = Z((Ai(i) - 1) + (1:length(mid))) ...
        + As(i) * mid;
end

[Bi, Bj, Bs] = find(INSTmodel.B);
for i = 1:length(Bi)   
    mid = [];
    for isol = 1:length(INSTmodel.Y{Bj(i)}{1})
        sxint = deval(INSTmodel.Y{Bj(i)}{1}(isol), t, INSTmodel.Y{Bj(i)}{2}{isol});
        mid = conv([ 1 - sum(mid); mid ], [ 1 - sum(sxint); sxint ]);
        mid = mid(2:end);
    end
    
    Z((Bi(i) - 1) + (1:length(mid))) = Z((Bi(i) - 1) + (1:length(mid))) ...
        + Bs(i) * mid;
end

Z = -INSTmodel.dC_diag ./ INSTmodel.C_diag .* Z;

for i = 1:length(Bi)   
    total_mid = 0;
    for isol = 1:size(INSTmodel.dY{Bj(i)}{1}, 1)
        mid = [];
        for jsol = 2:size(INSTmodel.dY{Bj(i)}{1}, 2)
            sxint = deval(INSTmodel.dY{Bj(i)}{1}(isol, jsol), t, INSTmodel.dY{Bj(i)}{2}{isol, jsol});
            mid = conv([ 1 - sum(mid); mid ], [ 1 - sum(sxint); sxint ]);
            mid = mid(2:end);
        end
        sxint = deval(INSTmodel.dY{Bj(i)}{1}(isol, 1), t, INSTmodel.dY{Bj(i)}{2}{isol, 1});
        mid = conv([ 1 - sum(mid); mid ], [ -sum(sxint); sxint ]);  % assumes first mid is a derivative
        mid = mid(2:end);
        
        total_mid = total_mid + mid;
    end
    
    Z((Bi(i) - 1) + (1:length(total_mid))) = Z((Bi(i) - 1) + (1:length(total_mid))) ...
        + Bs(i) * total_mid;
end

ddXdt = (INSTmodel.A * dX + Z) ./ INSTmodel.C_diag;
end