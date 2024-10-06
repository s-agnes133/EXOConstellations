import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from spherical_to_cartesian import calculate_distance_from_parallax, spherical_to_cartesian

# Assume 'calculate_distance_from_parallax' and 'spherical_to_cartesian' have been imported

# Function to calculate Cartesian coordinates for multiple stars
def plot_brightest_stars_3d(ra, dec, parallax):
    # Calculate distances from parallax (assuming parallax is in milliarcseconds)
    distances = np.array([calculate_distance_from_parallax(p) for p in parallax])

    # Convert RA, Dec, and Distance to Cartesian coordinates
    x, y, z = np.zeros(len(ra)), np.zeros(len(ra)), np.zeros(len(ra))
    for i in range(len(ra)):
        x[i], y[i], z[i] = spherical_to_cartesian(ra[i], dec[i], distances[i])
    
    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the stars as scatter points
    ax.scatter(x, y, z, c='blue', marker='o', s=50)

    # Set plot labels
    ax.set_xlabel('X (pc)')
    ax.set_ylabel('Y (pc)')
    ax.set_zlabel('Z (pc)')
    ax.set_title('3D Plot of 10 Brightest Stars')

    plt.show()

# Assuming 'brightest_ten' contains 'ra', 'dec', and 'parallax' for the ten brightest stars
# ra = brightest_ten['ra']
# dec = brightest_ten['dec']
# parallax = brightest_ten['parallax']

# Plot the stars in 3D
