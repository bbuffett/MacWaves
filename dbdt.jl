function  dbdt(r,x)

#   constructs sparse matrix Nsp for the time derivative equation
#
#   input
#   r        - radial grid (excludes boundaries)
#   x        - cos(theta) grid (excludes poles and equator)


   nr = length(r);
   nx = length(x);

   # initialize arrays
   ndim = nx * nr;
   row = Array{Int64}(undef,ndim);
   col = Array{Int64}(undef,ndim);
   val = Array{ComplexF64}(undef,ndim);

   #
   # time derivative of b_phi
   count = 0;
   for l = 1 : nx
       for k = 1 : nr

           count = count + 1;
           ipos = (l-1) * nr + k;
           row[count] = ipos;
           col[count] = ipos + ndim;
           val[count] = 1;
       end
   end

   nsp = sparse(row,col,val,2*ndim,2*ndim);

end
