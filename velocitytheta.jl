function velocitytheta(r,x,b,mode,d)

#   evaluates the theta component of velocity from the solution
#
#   input
#   r        - radial grid (excludes boundaries)
#   x        - cos(theta) grid (excludes poles and equator)
#   b        - radial magnetic field (mT)
#   mode     - solution vector (eigenfunction)
#   d        - eigenvalue  (frequency = d / 1im)


#   some constants
    omega = 0.7292e-4;       # rotation rate (1/s)
    sectoyear = 3.155e7;     # seconds in a year
    mu = 4 * pi * 1.0e-7;    # permeability
    eta = 0.8;            # magnetic diffusivity  (m^2/s)
    rho = 1.0e4;             # density
    L = 3.480e6;             # core radius (length scale)
    bscale = sqrt(rho * mu) * omega * L;

#   dimensionless quantities
    B = b * 1.0e-3 / bscale;
    Eta = eta / (omega * L^2);

#  allocate memory for velocity field
   nr = length(r);
   nx = length(x);
   ndim = nr * nx;
   dbdr = Array{ComplexF64}(undef,ndim);
   vthetacmb = Array{ComplexF64}(undef,nx)

#  radial magnetic field (dimensionless)
   br = magnetic_field(B,x);


#  evaluate radial gradient in bphi
   dr = r[2]-r[1];
   for l = 1 : nx
       for k = 2 : nr
           ipos = (l-1) * nr + k;
           dbdr[ipos] = (mode[ipos]-mode[ipos-1])/ dr;
       end
   end


#  extract CMB value
   for l = 1 : nx
       k = nr;
       ipos = (l-1)*nr + k;
       vthetacmb[l] = br[l]*dbdr[ipos]/(2*x[l]);
   end

   return vthetacmb

end
