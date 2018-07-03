function J = INSTfjac(t, y, INSTmodel)
%INSTFJAC calculates the jacobian of the ODEs solved in sqerror
%   J = INSTFJAC(T,Y,INSTMODEL) calculates the jacobian on the ODEs solved in
%   INST_sqerror.
%   T is time of the simulation
%   Y is the current value of the solution of the ODEs
%   INSTMODEL is a structure containing information on the ODEs tha describe
%   the current block. 

Cinv = spdiags(1 ./ INSTmodel.C_diag, 0, length(y), length(y));
J = Cinv * INSTmodel.A;
end