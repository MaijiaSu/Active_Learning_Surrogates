function MCMCchain = MCMCGenerator(Xstart,nsamples,pdf,proppdf,proprnd)      
    xi = Xstart;
    MCMCchain(1,:) = xi;
    for i = 1:nsamples-1

         xj = proprnd(xi);
         r = (pdf(xj)*proppdf(xj,xi))/(pdf(xi)*proppdf(xi,xj));

         if min([1,r]) > rand
             xi = xj;
         end

         MCMCchain(i+1,:) = xi;
    end
end