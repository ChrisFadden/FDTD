
function [field] = applyBC_E(field,grid)
%Magnetic Field Boundary Conditions
field.Ez(1) = field.Ez(2);

end
