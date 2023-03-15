import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
import numpy as np

proj = ccrs.LambertConformal(central_latitude = 25, 
                             central_longitude = 265, 
                             standard_parallels = (25, 25))

# Data and coordinates (from download link above)
with np.load('nam_218_20120414_1200_006.npz') as nam:
   dat = nam['dpc']
   lat = nam['lat']
   lon = nam['lon']
   
   
print (np.shape(lon))
print (np.shape(lat))
print (lon[0,:])
print (lat[0:-1:10,0])
print (lat[0:-1:10,1])
sys.exit()

ax = plt.axes(projection = proj)
ax.pcolormesh(lon, lat, dat, transform = ccrs.PlateCarree())
ax.add_feature(cf.NaturalEarthFeature(
               category='cultural',
               name='admin_1_states_provinces_lines',
               scale='50m',
               facecolor='none'))
ax.coastlines('50m')
ax.add_feature(cf.BORDERS)
plt.show()