function velocityphi(r,x,b,mode,d)

#   evaluates the velocity from the solution
#
#   input
#   r        - radial grid (excludes boundaries)
#   x        - cos(theta) grid (excludes poles and equator)
#   b        - radial magnetic field (mT)
#   eta      - dimensionless diffusivity
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
   dvdr = Array{ComplexF64}(undef,ndim);
   vphi = Array{ComplexF64}(undef,ndim);

#  radial magnetic field (dimensionless)
   br = magnetic_field(B,x);

#  retrieve required matrices
   dsp = diffusion(r,x,eta);

#  rate of change of bphi
   dbdt = d * mode[1:ndim];

#  account for diffusion
   dbdr2 = dsp[ndim+1:2*ndim,ndim+1:2*ndim] * mode[1:ndim];

#  combine to evaluate velocity gradient
   f = dbdt - dbdr2;
   for l = 1 : nx
       for k = 1 : nr
           ipos = (l-1) * nr + k;
           dvdr[ipos] = f[ipos] / br[l];
       end
   end

#  integrate for vphi
   dr = r[2] - r[1];
   for l = 1 : nx
       for k = 2 : nr
           ipos = (l-1) * nr + k;
           vphi[ipos] = vphi[ipos-1]+ 0.5*(dvdr[ipos] + dvdr[ipos-1])*dr;
       end
   end

#  extract CMB value
   vphicmb = Array{ComplexF64}(undef,nx)
   bphicmb = Array{ComplexF64}(undef,nx);

   for l = 1 : nx
       k = nr;
       ipos = (l-1)*nr + k;
#       println(" l ", l, " ipos ",ipos," vphi ",vphi[ipos]," ",vphi[ipos-1])
       vphicmb[l] = vphi[ipos];
       bphicmb[l] = mode[ipos];
#       println(" l ", vphicmb[l])
   end

   return vphicmb,bphicmb

end
