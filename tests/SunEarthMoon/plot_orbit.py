# Plot the orbit postion on the csv file
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

# Load the satellite positions from the CSV file
data = pd.read_csv('../../build/satellite_positions_sun_coords.csv')

# Plot the satellite positions in 3D
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot(data['X'], data['Y'], data['Z'], linestyle='-', color='b', label='Satellite Orbit', linewidth=0.2)
ax.set_xlabel('X (AU)')
ax.set_ylabel('Y (AU)')
ax.set_zlabel('Z (AU)')
ax.set_title('Satellite Orbit around the Moon')
ax.legend()

ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_zlim(-1.5, 1.5)

plt.grid()
plt.show()
