function  grids(nr,nx,rbottom)

#   constructs an equi-spaced grid in radius and x=cos(theta)
#   This particular version uses the shell thickness as a length
#   scale. This means that the radius at the top boundary is not r = 1;
#
#   nr       -  number of grid points in radius
#   nx       -  number of grid point in cos(theta)
#   rbottom   -  radius of lower layer  (dimensionless)
#
#   output
#   r        - radial grid (excludes boundaries)
#   x        - x grid (excludes poles and equator)
#   theta    - colatitude (radian)


   rcmb = 1.0;
   dr = (rcmb-rbottom) / nr;
   r = (rbottom+dr/2:dr:rcmb-dr/2);     # exclude boundaries

#  avoid grid point on equator
   if (rem(nx,2) ==1)
      dx = 2/(nx-1);
   else
       dx = 2/nx;
   end

   x = (-1 + dx/2 : dx : 1 - dx/2);
   theta = acos.(x);   # colatitude (in radian)

   return r, x, theta

end
