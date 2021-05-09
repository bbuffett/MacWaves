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
   br[1:nx] .= B;

#  dipole radial field
#  br = B * x;
#
#  increase rms field toward pole
#   br2 = B^2 * (0.48/0.65)^2 *(2.5*x.*x .+ 1);
#   br = sqrt.(br2);

   return br

end
