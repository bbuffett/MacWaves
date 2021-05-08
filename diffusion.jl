function diffusion(r,x,eta)

#   constructs sparse matrix Msp for "diffusion" or damping term
#
#   input
#   r        - radial grid (includes boundaries)
#   theta    - theta grid (excludes poles and equator)
#   eta      - dimensionless diffusivity


   nr = length(r);
   nx = length(x);
   dr = r[2]-r[1];
   dr2 = dr * dr;


#  initialize arrays
   ndim = nx * nr;      # number of grid points per variable
   nelements =  3 * ndim;
   row = Array{Int64}(undef,nelements);
   col = Array{Int64}(undef,nelements);
   val = Array{ComplexF64}(undef,nelements);


   #
   # fill elements of matrix (excluding boundaries)
   count = 0;
   for l = 1 : nx
       for k = 1 : nr

           ipos = ndim + (l-1) * nr +  k;

           # second derivative (central difference)
           if (k>1 && k < nr)
             count = count + 1;
             row[count] = ipos;
             col[count] = ipos;
             val[count] = -2 * eta /dr2;

             count = count + 1;
             row[count] = ipos;
             col[count] = ipos - 1;
             val[count] = 1 * eta /dr2;

             count = count + 1;
             row[count] = ipos;
             col[count] = ipos + 1;
             val[count] = 1 * eta /dr2;
           end

           if (k == 1)    # impose bc (uneven dr)

               count = count+1;
               row[count] = ipos;
               col[count] = ipos;
               val[count] = -5 * eta  / dr2;

               count = count + 1;
               row[count] = ipos;
               col[count] = ipos + 1;
               val[count] = 2 * eta /dr2;

               count = count + 1;
               row[count] = ipos;
               col[count] = ipos + 2;
               val[count] = (-1/5) * eta/ dr2;


           end

           if (k == nr)    # impose bc (uneven dr)

               count = count + 1;
               row[count] = ipos;
               col[count] = ipos;
               val[count] = -5 * eta / dr2;

               count = count + 1;
               row[count] = ipos;
               col[count] = ipos-1;
               val[count] = 2 * eta /dr2;

               count = count + 1;
               row[count] = ipos;
               col[count] = ipos - 2;
               val[count] = (-1/5) * eta / dr2;


           end


       end
   end


   dsp = sparse(row,col,val,2*ndim,2*ndim);


end
