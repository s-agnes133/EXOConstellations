import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# poprawić tę funkcję aby uwzględniała exoplanety
def calculate_distance_from_parallax(parallax_mas):
    """
    Calculate distance in parsecs from parallax in milliarcseconds (mas).
    
    Parameters:
    parallax_mas : float
        Parallax in milliarcseconds
    
    Returns:
    distance_pc : float
        Distance in parsecs
    """
    if parallax_mas is None or parallax_mas <= 0:
        return None  # Parallax should be positive for a valid distance
    
    distance_pc = 1000 / parallax_mas
    return distance_pc


# Function to convert RA, Dec, and Distance to Cartesian coordinates
def spherical_to_cartesian(ra, dec, distance):
    # Convert degrees to radians
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)
    
    # Cartesian coordinates
    x = distance * np.cos(dec_rad) * np.cos(ra_rad)
    y = distance * np.cos(dec_rad) * np.sin(ra_rad)
    z = distance * np.sin(dec_rad)
    
    return x, y, z

# # Example data for two stars (RA, Dec in degrees, Distance in parsecs)
# ra1, dec1, distance1 = 88.792939, 7.407064, 197  # Betelgeuse
# ra2, dec2, distance2 = 101.287155, -16.716116, 8.6  # Sirius

# # Convert to Cartesian coordinates
# x1, y1, z1 = spherical_to_cartesian(ra1, dec1, distance1)
# x2, y2, z2 = spherical_to_cartesian(ra2, dec2, distance2)

# # Create a 3D plot
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# # Plot the stars
# ax.scatter([x1, x2], [y1, y2], [z1, z2], color=['red', 'blue'], s=100, label=['Betelgeuse', 'Sirius'])

# # Labeling the stars
# ax.text(x1, y1, z1, 'Betelgeuse', color='red')
# ax.text(x2, y2, z2, 'Sirius', color='blue')

# # Labels and title
# ax.set_xlabel('X (pc)')
# ax.set_ylabel('Y (pc)')
# ax.set_zlabel('Z (pc)')
# ax.set_title('3D Position of Stars')

# plt.show()
