function macWave(nr,nx,b,Nmax,depth,period)

#   assembles the discretized equations and boundary conditions into
#   two matrices for a generalized eigenvalue problem. This version uses
#   the finite difference method on a wave equations for b_phi. The
#   solution does NOT include the diffusive magnetic boundary layer at the
#   based of the stratified layer as part of the solution. The equations
#   are expressed in non-dimensional form using diffusion-free scales
#
#   Scales
#   length      - radius
#   time        - 1/rotation rate  = 1 / omega
#   velocity    -  omega * L
#   mag. field  - sqrt(mu * rho) * omega * L
#
#   Inputs
#   nr      -  number of grid points in radius
#   nx      -  number of grid points in colatitude
#   b      -   radial (main) magnetic field (mT)
#   eta     -  magnetic diffusivity  (m^2/s)
#   Nmax   -   buoyancy frequency at CMB (Nmax/Omega)
#   depth   -  thickness of stratified layer (km)
#   period  -  target period for wave
#
#
#   Output
#   freq     -  eigenfrequency
#   q        -  quality factor
#   mode     -  eigenvector
#   r, x     -  r and x=cos(theta) coordinates
#
#
#   definition of eigenvalue problem
#   (A - d I) * v == 0
#

#   some constants
    omega = 0.7292e-4;       # rotation rate (1/s)
    sectoyear = 3.1557e7;     # seconds in a year
    mu = 4 * pi * 1.0e-7;    # permeability
    eta = 0.8;               # magnetic diffusivity  (m^2/s)
    rho = 1.0e4;             # density
    L = 3.480e6;             # core radius (length scale)
    bscale = sqrt(rho * mu) * omega * L;

#   dimensionless quantities
    B = b * 1.0e-3 / bscale;
    Eta = eta / (omega * L^2);

#     test for dipole and constant stratification
    eigenfreq1 = 0.5 * B * Nmax / 1.0;
    w1 = sqrt(1*2) * eigenfreq1
    w2 = sqrt(2*3) * eigenfreq1
    w3 = sqrt(3*4) * eigenfreq1
    w4 = sqrt(4*5) * eigenfreq1
    #println("w1 ",w1," w2 ", w2," w3 ",w3)


#   set up grid (shift grid to base of boundary layer)
    rbottom = (L - 1.0e3 * depth)/L;
    r,x,theta = grids(nr,nx,rbottom);

#   stratification
    iopt = 1;
    N2 = stratification(r,Nmax,rbottom,iopt);

#   magnetic field (constant in x or theta)
    br = magnetic_field(B,x);


#   set up linear system
     rsp = restoring(r,x,br,N2);
     dsp = diffusion(r,x,Eta);
     nsp = dbdt(r,x);
     asp = rsp + dsp + nsp;  # combine elements of A matrix


#  solve for eigenvalues ()
    freq = 2 * pi / (period * sectoyear * omega)
    shift = freq * 1im;
    d,v = eigs(asp;which=:LM,nev=8,sigma=shift,maxiter=600);


    return r,x,d,v


end
