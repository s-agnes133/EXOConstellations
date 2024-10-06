import pandas as pd
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import plotly.graph_objs as go
from flask import Flask, render_template, request
import logging

app = Flask(__name__)

# Ustawienia logowania
logging.basicConfig(level=logging.DEBUG)

def get_exoplanet_data():
    url = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_name,ra,dec,pl_bmasse+from+ps&format=csv'
    exoplanets = pd.read_csv(url)
    return exoplanets.head(30)  # Tylko 30 najpopularniejszych egzoplanet

def get_gaia_data(ra, dec, radius=5):
    coord = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')
    Gaia.ROW_LIMIT = 50000  # Limit gwiazd
    job = Gaia.cone_search_async(coordinate=coord, radius=(radius * u.degree))
    results = job.get_results()

    logging.debug(f"Data from Gaia: {results}")
    return pd.DataFrame({
        'ra': results['ra'],
        'dec': results['dec'],
        'parallax': results['parallax'],
        'phot_g_mean_mag': results['phot_g_mean_mag']
    })

def transform_coordinates(gaia_data, exoplanet_coord):
    valid_data = gaia_data[gaia_data['parallax'] > 0].copy()
    if valid_data.empty:
        raise ValueError("No valid parallax data available for transformation.")

    # Tworzenie obiektu współrzędnych nieba dla gwiazd
    sky_coords = SkyCoord(
        ra=valid_data['ra'].values * u.deg,
        dec=valid_data['dec'].values * u.deg,
        distance=(1 / valid_data['parallax'].values) * u.pc,  # Konwersja paralaksy na odległość w parsekach
        frame='icrs'
    )

    logging.debug(f"Sky coordinates: {sky_coords}")

    # Obliczamy różnice w odległości między egzoplanetą a każdą gwiazdą
    transformed_coords = sky_coords.transform_to(SkyCoord(ra=exoplanet_coord.ra, dec=exoplanet_coord.dec, frame='icrs'))

    return pd.DataFrame({
        'ra_transformed': transformed_coords.ra.deg,
        'dec_transformed': transformed_coords.dec.deg,
        'distance': transformed_coords.distance.pc  # Odległości w parsekach
    })

def scale_and_center_coordinates(data, scaling_factor=1e3):
    """Przekształca RA, DEC i odległość na współrzędne kartezjańskie oraz przesuwa egzoplanetę do centrum."""

    # Konwersja RA, DEC i odległości na współrzędne kartezjańskie (x, y, z)
    ra_rad = np.deg2rad(data['ra_transformed'])  # RA w radianach
    dec_rad = np.deg2rad(data['dec_transformed'])  # DEC w radianach
    distance = data['distance'].values  # Odległość w parsekach

    # Przekształcenie współrzędnych sferycznych na kartezjańskie
    x = distance * np.cos(dec_rad) * np.cos(ra_rad)
    y = distance * np.cos(dec_rad) * np.sin(ra_rad)
    z = distance * np.sin(dec_rad)

    # Przesunięcie układu tak, aby egzoplaneta była w centrum
    x_centered = x - np.mean(x)
    y_centered = y - np.mean(y)
    z_centered = z - np.mean(z)

    # Skala, aby wartości nie były zbyt duże
    x_centered /= scaling_factor
    y_centered /= scaling_factor
    z_centered /= scaling_factor

    return pd.DataFrame({
        'x_centered': x_centered,
        'y_centered': y_centered,
        'z_centered': z_centered
    })


def create_visualization(transformed_data,exoplanet_names, scaling_factor=1e3):
    """Tworzy wizualizację 3D gwiazd wokół egzoplanety z identycznymi zakresami osi."""

    # Skalowanie i centrowanie współrzędnych
    transformed_data = scale_and_center_coordinates(transformed_data, scaling_factor)

    # Punkt dla egzoplanety w centrum (0, 0, 0)
    exoplanet_trace = go.Scatter3d(
        x=[0], y=[0], z=[0],  # Egzoplaneta w punkcie (0, 0, 0)
        mode='markers',
        marker=dict(
            size=5,  # Zmniejszenie wielkości punktu egzoplanety
            color='red',
            opacity=1.0
        ),
        name=str(exoplanet_names[0])
    )

    # Punkty dla gwiazd wokół egzoplanety
    stars_trace = go.Scatter3d(
        x=transformed_data['x_centered'],
        y=transformed_data['y_centered'],
        z=transformed_data['z_centered'],
        mode='markers',
        marker=dict(
            size=2,
            color='blue',  # Kolor gwiazd
            opacity=0.8
        ),
        name='Objects'
    )

    # Ustawienie identycznych zakresów dla osi X, Y, Z
    max_range = max(
        transformed_data['x_centered'].max() - transformed_data['x_centered'].min(),
        transformed_data['y_centered'].max() - transformed_data['y_centered'].min(),
        transformed_data['z_centered'].max() - transformed_data['z_centered'].min()
    )

    # Zakresy na osiach z wartościami ujemnymi, aby egzoplaneta była w centrum
    x_range = [-max_range / 2, max_range / 2]
    y_range = [-max_range / 2, max_range / 2]
    z_range = [-max_range / 2, max_range / 2]

    # Ustawienia layoutu dla 3D z równomiernymi podziałkami na osiach
    layout = go.Layout(
        title='View in 3D space around an exoplanet',
        scene=dict(
            xaxis=dict(
                title='X (pc)',  # Dodanie opisu osi
                backgroundcolor='white',
                gridcolor='gray',
                showbackground=True,
                range=x_range,  # Zastosowanie identycznego zakresu
                zeroline=True,
                zerolinecolor='gray'
            ),
            yaxis=dict(
                title='Y (pc)',  # Dodanie opisu osi
                backgroundcolor='white',
                gridcolor='gray',
                showbackground=True,
                range=y_range,  # Zastosowanie identycznego zakresu
                zeroline=True,
                zerolinecolor='gray'
            ),
            zaxis=dict(
                title='Z (pc)',  # Dodanie opisu osi
                backgroundcolor='white',
                gridcolor='gray',
                showbackground=True,
                range=z_range,  # Zastosowanie identycznego zakresu
                zeroline=True,
                zerolinecolor='gray'
            ),
            aspectmode='cube'  # Zachowanie proporcji jako sześcian
        ),
        paper_bgcolor='white',  # Białe tło całego wykresu
        plot_bgcolor='white'  # Białe tło samego obszaru wykresu
    )

    # Tworzenie wykresu z egzoplanetą w centrum i gwiazdami dookoła
    fig = go.Figure(data=[exoplanet_trace, stars_trace], layout=layout)
    return fig.to_html(full_html=False)


@app.route('/')
def index():
    exoplanets = get_exoplanet_data()
    return render_template('index.html', exoplanets=exoplanets)

@app.route('/visualize', methods=['POST'])
def visualize():
    logging.debug("Wizualizacja rozpoczęta.")
    logging.debug(f"Request form data: {request.form}")

    exoplanet_names = request.form.getlist('exoplanet')
    logging.debug(f"Wybrane egzoplanety: {exoplanet_names}")

    if not exoplanet_names:
        return "Proszę wybrać przynajmniej jedną egzoplanetę.", 400

    visualizations = []

    for exoplanet_name in exoplanet_names:
        # Pobieranie danych wybranej egzoplanety
        selected_exoplanet = get_exoplanet_data().set_index('pl_name').loc[exoplanet_name]

        ra = selected_exoplanet['ra']
        dec = selected_exoplanet['dec']

        # Przekształcanie współrzędnych egzoplanety
        exoplanet_coord = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')

        # Pobieranie danych Gaia dla danej lokalizacji
        gaia_data = get_gaia_data(ra, dec)

        # Sprawdzenie, czy pobrane dane Gaia są puste
        if gaia_data.empty:
            logging.error("Brak danych Gaia dla egzoplanety.")
            return "Brak danych Gaia dla wybranej egzoplanety.", 400

        # Transformacja współrzędnych do układu odniesienia egzoplanety
        transformed_data = transform_coordinates(gaia_data, exoplanet_coord)

        # Tworzenie wizualizacji z egzoplanetą w centrum
        visualization_html = create_visualization(transformed_data,exoplanet_names)
        visualizations.append(visualization_html)

    return render_template('visualization.html', visualizations=visualizations)

if __name__ == '__main__':
    app.run(debug=True)