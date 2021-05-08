function inertia(r,x)

#   constructs sparse matrix Msp for "inertial" term
#
#   input
#   r        - radial grid (includes boundaries)
#   x        - cos(theta) grid (excludes poles and equator)
#
#  output
#  msp       - sparse m matrix


   nr = length(r);
   nx = length(x);
   c =  1.0 ;  #complex(0,1);

#  initialize arrays
   ndim = nx * nr;
   nbc = 2 * nx;
   nelements = 2 * ndim;
   row = Array{Int64}(undef,nelements);
   col = Array{Int64}(undef,nelements);
   val = Array{ComplexF64}(undef,nelements);
  # val = Array{Float64}(undef,nelements);

#   row = zeros(nelements);
#   col = zeros(nelements);
#   val = zeros(nelements);

#
#  time derivative of x = b_phi
   count = 0;
   for l = 1 : nx
       for k = 1 : nr

           count = count + 1;
           ipos = (l-1) * nr + k;
           row[count] = ipos;
           col[count] = ipos;
           val[count] = c;
       end
   end

#  time derivative of y = d_t phi
   for l = 1 : nx
       for k = 1 : nr

           count = count + 1;
           ipos = ndim + (l-1) * nr + k;
           row[count] = ipos;
           col[count] = ipos;
           val[count] = c;
       end
   end

   msp = sparse(row,col,val,2*ndim+nbc,2*ndim+nbc);

end
