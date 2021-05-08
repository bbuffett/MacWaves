function stratification(r,Nmax,rbottom,iopt)

# defines stratification on radial grid.

# r       - radial grid points
# Nmax    - stratification at CMB (r=1)
# Rbottom - radius of base of stratified layer
# iopt =  - constant stratification (=0); linear stratification (=1)

#   constructs density stratification
   n = length(r);
   N2 = zeros(n);



   N2_peak = Nmax * Nmax;


#  stratified layer
   if (iopt == 1)
       N2[1:n] .= (r[1:n] .- r[1]) * N2_peak / (r[n] - r[1]);
   end

   if (iopt == 0 )
       N2[1:n] .= N2_peak;
   end


   return N2

end
