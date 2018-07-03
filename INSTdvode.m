function ddXdt = INSTdvode(t, dX, INSTmodel)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Z = zeros(size(dX));
[Bi, Bj, Bs] = find(INSTmodel.B);
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

dA = INSTmodel.dA(:, cumsum([1; INSTmodel.emu_natoms(1:(end - 1))]));
[dAi, dAj, dAs] = find(dA);
for i = 1:length(dAi)   
    mid = [];
    for isol = 1:length(INSTmodel.X{dAj(i)}{1})
        sxint = deval(INSTmodel.X{dAj(i)}{1}(isol), t, INSTmodel.X{dAj(i)}{2}{isol});
        mid = conv([ 1 - sum(mid); mid ], [ 1 - sum(sxint); sxint ]);
        mid = mid(2:end);
    end
    
    Z((dAi(i) - 1) + (1:length(mid))) = Z((dAi(i) - 1) + (1:length(mid))) ...
        + dAs(i) * mid;
end

[dBi, dBj, dBs] = find(INSTmodel.dB);
for i = 1:length(dBi)   
    mid = [];
    for isol = 1:length(INSTmodel.Y{dBj(i)}{1})
        sxint = deval(INSTmodel.Y{dBj(i)}{1}(isol), t, INSTmodel.Y{dBj(i)}{2}{isol});
        mid = conv([ 1 - sum(mid); mid ], [ 1 - sum(sxint); sxint ]);
        mid = mid(2:end);
    end
    
    Z((dBi(i) - 1) + (1:length(mid))) = Z((dBi(i) - 1) + (1:length(mid))) ...
        + dBs(i) * mid;
end

ddXdt = (INSTmodel.A * dX + Z) ./ INSTmodel.C_diag;
end