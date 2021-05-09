function  magnetic_field(B,x)

# evaluates the radial magnetic field as a function of x = cos(theta);
#
# input
# B   - amplitude of magnetic field (dimensionless)
# x   - cos(theta)
#
# output
# br     - radial field as a function of x

#  size of x grid
   nx = length(x);
   br = zeros(nx);

#  constant radial field
#   br[1:nx] .= B;

#  dipole radial field
#  br = B * x;
#
#  increase rms field toward pole
   Gamma = 2.5;
   B02 = B^2 / (1.0  + Gamma/3.0);
   br2 = B02 *(Gamma*x.*x .+ 1);
   br = sqrt.(br2);

   return br

end
