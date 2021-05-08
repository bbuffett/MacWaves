function flow(t,t0,lat,long,x,vphi,d)

# computes the velocity and acceleration field on a lat,long
# grid at time t using a mode from macWave()
#
# input
#
# t    - time (years)
# t0   - origin time (years)
# colat - latitude (degree)
# long - longitude (degree)
# x    - meridional coordinate (x = cos(theta))
# vphi - solution for vphi at CMB (vs. x)
# d    - eigenvalue ( i x frequency )

#
# output
#
# vx   - phi component of velocity at CMB
# ax   - phi components of acceleration

# some constants
sectoyear = 365.25 * 24 * 60 * 60;
tsec = (t-t0) * sectoyear;
omega = 0.7292e-4;

# add solution at poles x = -1, 1
n = length(x)
vphi_ends = Array{ComplexF64}(undef,n+2);
x_ends = Array{Float64}(undef,n+2);
x_ends[1] = -1.0;
x_ends[2:n+1] = x;
x_ends[n+2] = 1.0;
vphi_ends[1] = 0.0;
vphi_ends[2:n+1] = vphi;
vphi_ends[n+2] = 0.0;

# allocate memory for velocity and acceleration on colat grid
vr = Array{Float64,1}(undef,length(lat));
vi = Array{Float64,1}(undef,length(lat));
v = Array{Float64,2}(undef,length(lat),length(long));
a = Array{Float64,2}(undef,length(lat),length(long));

# dimensional eigenvalue
dd = d * omega;

# evaluate latitude of wave solution
theta = -acos.(x_ends).+pi/2;

#  interpolate wave on to lat,long grid
nodes = (theta,);
itp = Interpolations.interpolate(nodes,real(vphi_ends),Gridded(Linear()));
vr = itp(lat);
itp = Interpolations.interpolate(nodes,imag(vphi_ends),Gridded(Linear()));
vi = itp(lat);

# copy real part onto grid
for j = 1 : length(lat)
    for k = 1 : length(long)

       phase = exp(dd * tsec);
       v[j,k] = real( (vr[j] + vi[j]*1im)*phase);
       a[j,k] = real( (vr[j] + vi[j]*1im)*dd*phase)*sectoyear;

    end
end

# normalize
vmax = maximum(v);
v = v / vmax;
a = a / vmax;

return v,a

end
