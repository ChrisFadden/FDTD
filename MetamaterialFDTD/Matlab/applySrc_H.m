function [field] = applySrc_H(field,grid,src,constants,t)

%TF/SF Source Implementation
offset = (grid.TFSF_end - grid.TFSF_start)*grid.dt;

%Gaussian Source
srcStart = exp(-((t*grid.dt - src.t0)/src.tw)^2);
srcEnd = exp(-((t*grid.dt - src.t0 - offset)/src.tw)^2);

%Harmonic Source (need no u(t - tfsf) since the "error" is just in
%                 scattered field, it won't reflect off a structure                
%srcStart = sin(2*pi*t*grid.dt*src.freq);
%srcEnd = sin(2*pi*src.freq*(t*grid.dt - offset));

%Field Updates
field.Hy(grid.TFSF_start) = field.Hy(grid.TFSF_start) + ...
                            sqrt(constants.eps0 / constants.mu0)*...
                            srcStart;
field.Hy(grid.TFSF_end) = field.Hy(grid.TFSF_end) - ...
                             sqrt(constants.eps0 / constants.mu0)*...
                             srcEnd;    
  

end

