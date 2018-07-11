#-------------------------------------------------------------------------
#					IDV_GRIDUTIL.PY
#--------------------------------------------------------------------------
# The following functions were adapted from the GridUtil.java code from IDV software
# by Leah Campbell, 6 October 2016
# Primarily made for smoothing functions, although there may be more in GridUtil.java

#------------------------------------------------------------------------------------------------
'''
1. SMOOTHER FUNCTION 
 * Apply a weigthed smoothing function to the grid.  The smoothing types are:
 * <p>
 * SMOOTH_CRESSMAN: the smoothed value is given by a weighted average of values
 * at surrounding grid points.  The weighting function is the Cressman weighting
 * function:
 * <pre>
 *         w = ( D**2 - d**2 ) / ( D**2 + d**2 )
 * </pre>
 * In the above, d is the distance (in grid increments) of the neighboring point
 * to the smoothing point, and D is the radius of influence [in grid increments]
 * <p>
 * SMOOTH_CIRCULAR: the weighting function is the circular apperture
 * diffraction function (following a suggestion of Barnes et al. 1996):
 * <pre>
 *          w = bessel(3.8317*d/D)/(3.8317*d/D)
 * </pre>
 * <p>
 * SMOOTH_RECTANGULAR: the weighting function is the product of the rectangular
 * apperture diffraction function in the x and y directions (the function used
 * in Barnes et al. 1996):
 * <pre>
 *          w = [sin(pi*x/D)/(pi*x/D)]*[sin(pi*y/D)/(pi*y/D)]
 * </pre>
 * Apply a rectangular aperature smoothing to the grid points.  The weighting 
 function is the product of the rectangular aperature diffraction function
 in the x and y directions.  D is the radius of influence in grid 
 increments, increasing D increases the smoothing. (default D=2)
 
 * Adapted from smooth.f written by Mark Stoelinga in his RIP package
 * 
'''

# function options:
# orig_array = original 2d array of data (haven't tried with 3d....)
# smoothing_scheme = 'cressman', 'rect', or 'circular'
# radius = radius of window in grid units, default is 2
#

#el_file = '/uufs/chpc.utah.edu/common/home/u1013082/elevation_data/60second_westernUS.nc'
#fh = Dataset(el_file, mode='r')
#elevation_p = fh.variables['Band1'][:]
#lat_p = fh.variables['lat'][:]
#lon_p = fh.variables['lon'][:]
#
#xlen = len(elevation_p[:,0])
#ylen = len(elevation_p[0,:])
#
#lat_netcdf_p = zeros((xlen, ylen))
#long_netcdf_p = zeros((xlen,ylen))
#for i in range(xlen):
#    long_netcdf_p[i,:] = lon_p
#for j in range(ylen):
#    lat_netcdf_p[:,j] = lat_p
    
 


def smoother(orig_array,smoothing_scheme,radius):
	import numpy as np
	import pdb

	# return original array if radius is 0
	if radius == 0:
		return orig_array
	
	maxwts = 100 # maximum number of weights
	nfp = min(maxwts, 2*radius) # can't have a radius bigger than 50
	npsq = float(radius * radius) # used in Cressman scheme
	beszero = 3.8317 # used in Circular scheme
	orig_array = np.array(orig_array)
	pslab = orig_array.flatten() # flattened original array
	newvalues = np.empty(np.shape(pslab))	# empty array to fill with smoothed values
	fprt = np.empty([nfp,nfp])	# smoother array, changes depending on smoother option
	
	# create smoother array
	if smoothing_scheme == 'rect': # for rectangular aperture
		for i in range(0,nfp):
			for j in range(0, nfp):
				if j == radius:
					xfac = 1.0
				else:
					xdist = np.pi / radius * (j - radius)
					xfac = np.sin(xdist) / xdist
				if i == radius:
					yfac = 1.0
				else:
					ydist = np.pi / radius * (i - radius)
					yfac  = np.sin(ydist) / ydist
				fprt[j,i] = xfac * yfac
				
	elif smoothing_scheme == 'cressman': # For Cressman scheme
		for i in range(0,nfp):
			for j in range(0, nfp):
				
				distsq = np.power((i - radius), 2) + np.power((j-radius), 2)
				fprt[j,i] = max((npsq - distsq) / (npsq + distsq), 0.0)
			
	elif smoothing_scheme == 'circular': # For circular diffraction scheme
		for i in range(0,nfp):
			for j in range(0, nfp):
				dist = (beszero / radius * np.sqrt( np.power((i - radius), 2) + np.power((j - radius), 2)))
				if (i == radius) and (j == radius):
					fprt[j,i] = 0.5
				else:
					fprt[j,i] = max(0.0, bessel(dist) / dist)

	#now do the work of smoothing
	niy = np.shape(orig_array)[0]
	njx = np.shape(orig_array)[1]
	for i in range(0, niy):
		for j in range(0,njx):
			#pdb.set_trace()
			index = j + i * njx
			if not np.isnan(pslab[index]):
				tot = 0.0
				totwt = 0.0
				iss = max(0, i - radius)
				ie = min(niy - 1, i + radius)
				jss = max(0, j - radius)
				je = min(njx - 1, j + radius)
				for ireg in range(iss, ie):
					ifp = ireg - i + radius
					for jreg in range(jss, je):
						jfp = jreg - j + radius
						psindex = ireg * njx + jreg
						if not np.isnan(pslab[psindex]):
							totwt = totwt + fprt[jfp,ifp]
							tot = tot + fprt[jfp,ifp] * pslab[psindex]
				try:
					newvalues[index] = tot / totwt
				except:
					newvalues[index] = np.nan
			else:
				newvalues[index] = np.nan

	return np.reshape(newvalues,[niy,njx])

#new_array = smoother(elevation_p, 'rect', 3)

#------------------------------------------------------------------------------------------------
'''
2. BESSEL FUNCTION 
used in Circular smoothing scheme
Copied from RIP
'''
# x = the value
# returns the function

def bessel(x):
	import numpy as np
	
	rint = 0.0
	for i in range(0,1000):
		u = i * 0.001 - 0.0005
		rint = rint + np.sqrt(1.0 - u * u) * np.cos(x * u) * 0.001
	
	return 2.0 * x * rint / (4.0 * np.arctan(10))

#------------------------------------------------------------------------------------------------

'''
3. HORIZONTAL DIVERGENCE
Copied from RIP
'''
'''
def div(V):
  """ Horizontal Divergence 
  <div class=jython>
      DIV ( V ) = DDX ( u ) + DDY ( v ) 
  </div>
  """
  return add(ddx(ur(V)),ddy(vr(V)))
def ddx(S):
  """ Take the derivative with respect to the domain's X coordinate 
  """
  return GridMath.ddx(S);
  
def ur(V):
  """ Grid relative u component 
  """
  return DerivedGridFactory.getUComponent(V)
  
def vr(V):
  """ Grid relative v component 
  """
  return DerivedGridFactory.getVComponent(V)
'''
