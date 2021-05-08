function  restoring(r,x,br,N2)

#   constructs part of sparse matrix wave operator Lsp for the
#   restoring force. It involves the L_x^2 term in the wave eqn.
#
#
#   input
#   r       - radial grid (includes boundaries)
#   x       - latitude grid x = cos(theta) (excludes poles and equator)
#   br      - dimensionless radial field
#   N2      - dimensionless stratification (N/Omega)^2
#


   nr = length(r);
   nx = length(x);
   dx = x[2]-x[1];
   dx2 = dx * dx;
   rcmb = 1.0;

   # some useful arrays
   s = sqrt.(1 .- x.*x);
   t = s./x;
   b = 0.25 * br / (rcmb^2);   # less useful now that L = Rcmb


   # initialize arrays
   ndim = nx * nr;
   nelements = 3 * ndim;
   row = Array{Int64}(undef,nelements);
   col = Array{Int64}(undef,nelements);
   val = Array{ComplexF64}(undef,nelements);

#println("nr ",nr," nx ",nx)
#
#  fill matrix (including boundaries)
   count = 0;
   for l = 1 : nx
       for k = 1 : nr

           ipos = (l-1) * nr +  k;

           # second derivative in x (central difference)
           if (l>1 && l < nx)


             count = count + 1;
             row[count] = ndim + ipos;
             col[count] = ipos-nr;
             val[count] = t[l] * b[l] * N2[k] * br[l-1] * t[l-1] /dx2;

             count = count + 1;
             row[count] = ndim + ipos;
             col[count] = ipos;
             val[count] = -2*t[l] * b[l] * N2[k] * br[l] * t[l] /dx2;

             count = count + 1;
             row[count] = ndim + ipos;
             col[count] = ipos + nr;
             val[count] = t[l] * b[l] * N2[k] * br[l+1] * t[l+1] /dx2;

           end

           if (l == 1)   # center difference with pole condition; unequal dx)

               count = count+1;
               row[count] = ndim + ipos;
               col[count] = ipos;
               val[count] = -5*t[l]*b[l]*N2[k] * br[l] * t[l] /dx2;

               count = count + 1;
               row[count] = ndim + ipos;
               col[count] = ipos + nr;
               val[count] = (2)*t[l]*b[l]*N2[k] * br[l+1] * t[l+1] /dx2;

               count = count + 1;
               row[count] = ndim + ipos;
               col[count] = ipos + 2*nr;
               val[count] = -(1/5)*t[l]*b[l]*N2[k]* br[l+2] * t[l+2]/dx2;

           end

           if (l == nx)   # backward derivative

               count = count + 1;
               row[count] = ndim + ipos;
               col[count] = ipos;
               val[count] = -5*t[l] * b[l] * N2[k] * br[l] * t[l] /dx2;

               count = count + 1;
               row[count] = ndim + ipos;
               col[count] = ipos-nr;
               val[count] = (2)*t[l]*b[l]*N2[k] * br[l-1] * t[l-1] /dx2;

               count = count + 1;
               row[count] = ndim + ipos;
               col[count] = ipos - 2*nr;
               val[count] = -(1/5)*t[l]*b[l]*N2[k] * br[l-2]* t[l-2]/dx2;

           end


       end
   end

#  add boundary conditions on top and bottom radii
   rsp = sparse(row,col,val,2*ndim,2*ndim);






end
