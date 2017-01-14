function [field] = applyBC_H(field,grid)
%Magnetic Field Boundary Conditions

field.Hy(grid.sizeX) = field.Hy(grid.sizeX-1);

end

