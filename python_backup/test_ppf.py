import scipy.stats as stats


rq = 1.0
rp = 3.5
        
# ---- the number of (interior) knots in the cubic spline will
#      be set to be no greater than 9.0, and for a sample size
#      of 100 will be 3.
        
nz = 205
nknots = min([9,nz//30])
print ('nknots = ', nknots)
query_these_indices = []
cdf_at_indices = []
        
for iknot in range (1,nknots+1):
    rknot = float(iknot)/(nknots+1)
    print ('rknot, rp, rq = ', rknot, rp, rq)
    xloc = stats.beta.ppf(rknot, rp, rq)
    c = stats.beta.cdf(xloc, rp, rq)
    print ('c = ', c)
    iloc = int(nz*xloc)
    print ('nknots, iknot, rknot, xloc, iloc, nz = ', \
        nknots, iknot, rknot, xloc, iloc, nz)
    query_these_indices.append(iloc)
    c = (1./(2.*nz)) + float(iloc)/float(nz)
    cdf_at_indices.append(c)
        
print ('nz = ', nz)    
print ('nknots = ', nknots)
print ('query_these_indices = ', query_these_indices)
print ('cdf_at_indices = ', cdf_at_indices)
        
