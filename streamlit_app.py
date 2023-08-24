import streamlit as st
from streamlit_folium import folium_static
import folium
import rasterio
from pysheds.grid import Grid
import numpy as np
import matplotlib.pyplot as plt
import math



def haversine_distance(coord1, coord2):
    R = 6371.0  # Earth radius in kilometers

    lat1, lon1 = coord1
    lat2, lon2 = coord2

    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)

    a = math.sin(dlat / 2) ** 2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    distance = R * c
    return distance


def calculate_area(coordinates):
    total_area = 0
    num_points = len(coordinates)

    for i in range(num_points - 1):
        total_area += haversine_distance(coordinates[i], coordinates[i + 1])
        # Area in square kilometers
    return total_area

st.set_page_config(
    page_title="Hydrology App",
    page_icon="🧊",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        'About': "Copyright of Anurak Puengrostham"
    }
)

digital_terrain = "data/mergeDEM.tif"
streams = 'data/export.tif'
@st.cache_data()
def grid_from_raster(file):
    grid = Grid.from_raster(file)
    return grid

@st.cache_resource()
def open_rasterio(file):
    dem = rasterio.open(file)
    return dem

grid = grid_from_raster(digital_terrain)
grid_array = grid.read_raster(digital_terrain)

digital_terrain_src = open_rasterio(digital_terrain)
digital_terrain_array = digital_terrain_src.read()

stream_raster = rasterio.open(streams)
array = stream_raster.read()
bounds = digital_terrain_src.bounds
bbox = [(bounds.bottom, bounds.left), (bounds.top, bounds.right)]

st.title("Watershed Delineation")
st.text("by Anurak Puengrostham (Copyright)")
st.text(" ")
st.text("Source of data:-")
st.text("1. Digital Elevation from Earth Explorer USGS")
st.text("2. DEM is processed to remove pits, depression, and flatness")
st.text("3. Area of interest is limited to Project Area only!")

st.header('1. Input Pour Point:')
st.text(f"Input Pour-point - from Lat:{bounds.bottom} to {bounds.top}\n"
        f"Long: from {bounds.left} to {bounds.right}")
lat = st.number_input('Insert Latitude', min_value= bounds.bottom, max_value= bounds.top, value= 19.889)
st.write('The current number is ', lat)
long = st.number_input('Insert Longtitude', min_value=bounds.left, max_value=bounds.right, value= 99.878)
st.write('The current number is ', long)

# center on the map
n = folium.Map(location=[lat,long], zoom_start=14)

tooltip = f"Pour Point: {lat},{long}"
folium.Marker(
    [lat, long], popup="Double Track Project", tooltip=tooltip
).add_to(n)
folium.LayerControl().add_to(n)
folium_static(n)

terrain_button = st.button("Digital Terrain MAP")
stream_button = st.button("Streams Line")
catchment_button = st.button("Catchment")
if terrain_button:
    st.header('Digital Terrain Map')
    # center on the map
    n = folium.Map(location=[lat,long], zoom_start=10)

    tooltip = f"Pour Point: {lat},{long}"
    folium.Marker(
        [lat, long], popup="Double Track Project", tooltip=tooltip
    ).add_to(n)

    terrain_img = folium.raster_layers.ImageOverlay(
        name="DEM",
        image=np.moveaxis(digital_terrain_array, 0, -1),
        bounds=bbox,
        opacity=0.7,
        interactive=True,
        cross_origin=False,
        zindex=1,
    )
    # folium.Popup("Message").add_to(img)
    terrain_img.add_to(n)
    folium.LayerControl().add_to(n)
    folium_static(n)
else: pass

if stream_button:
    st.header('Stream Lines Map')
    # center on the map
    m = folium.Map(location=[lat,long], zoom_start=10)

    # add marker for Liberty Bell
    tooltip = f"lat:long = [{lat}:{long}]"
    folium.Marker(
        [lat, long], popup=f"[{lat}:{long}]", tooltip=tooltip
    ).add_to(m)

    img = folium.raster_layers.ImageOverlay(
        name="Stream",
        image=np.moveaxis(array, 0, -1),
        bounds=bbox,
        opacity=0.5,
        interactive=True,
        cross_origin=False,
        zindex=1,
    )

    # folium.Popup("Message").add_to(img)
    img.add_to(m)
    folium.LayerControl().add_to(m)
    folium_static(m)
else: pass

if catchment_button:
    st.header('Catchment Area')
    # Condition DEM
    # ----------------------
    # Fill pits in DEM
    pit_filled_dem = grid.fill_pits(grid_array)

    # Fill depressions in DEM
    flooded_dem = grid.fill_depressions(pit_filled_dem)

    # Resolve flats in DEM
    inflated_dem = grid.resolve_flats(flooded_dem)

    # Determine D8 flow directions from DEM
    # ----------------------
    # Specify directional mapping
    dirmap = (64, 128, 1, 2, 4, 8, 16, 32)
    # Compute flow directions
    # -------------------------------------
    fdir = grid.flowdir(inflated_dem, dirmap=dirmap)
    # -------------------------------------
    # Calculate flow accumulation
    # --------------------------
    acc = grid.accumulation(fdir, dirmap=dirmap)
    # Specify pour point
    print(f'lat, long :{long},{lat}')

    # Snap pour point to high accumulation cell
    accum = 500 # Initial set variable
    x_snap, y_snap = grid.snap_to_mask(acc > accum, (long, lat))

    # Delineate the catchment
    catch = grid.catchment(x=x_snap, y=y_snap, fdir=fdir, dirmap=dirmap,
                           xytype='coordinate')

    # Crop and plot the catchment
    # ---------------------------
    # Clip the bounding box to the catchment
    grid.clip_to(catch)
    clipped_catch = grid.view(catch)
    lat1,lat2, long1, long2 = clipped_catch.extent
    coor = [
        (lat1,long1),
        (lat1,long2),
        (lat2,long2),
        (lat2,long1),
        (lat1,long1)
    ]
    area = calculate_area(coor)
    total_cells = clipped_catch.shape[0] * clipped_catch.shape[1]
    catchment_percent = clipped_catch.sum()/total_cells
    st.write("Optimum Pour Point:", f'{x_snap:.3f}', f'{y_snap:.3f}')
    st.write("Catchment Area = ", f'{(area * catchment_percent):.2f}', "Square Kilometers")
    st.write("Latidude  : distance :  ", f'{haversine_distance((lat1, long1), (lat2, long1)):.2f}', "  km",
             "Cell :", clipped_catch.shape[0])
    st.write("Longitude : distance :  ", f'{haversine_distance((lat1, long1), (lat1, long2)):.2f}', "  km",
             "Cell :", clipped_catch.shape[1])

    # Calculate distance to outlet from each cell
    # -------------------------------------------
    dist = grid.distance_to_outlet(x=x_snap, y=y_snap, fdir=fdir, dirmap=dirmap, xytype='coordinate')

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_alpha(0)
    plt.grid('on', zorder=0)
    im = ax.imshow(dist, extent=grid.extent, zorder=2,
                   cmap='cubehelix_r')
    plt.colorbar(im, ax=ax, label='Distance to outlet (cells)')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Flow Distance', size=14)
    st.pyplot(fig)

    # Catchment imposed on Terrain
    catchment_on_terrain = np.where(np.invert(catch), digital_terrain_array,0)

    # center on the map
    n = folium.Map(location=[lat, long], zoom_start=10)

    tooltip = f"Pour Point: {lat},{long}"
    folium.Marker(
        [lat, long], popup="Double Track Project", tooltip=tooltip
    ).add_to(n)

    terrain_img = folium.raster_layers.ImageOverlay(
        name="DEM",
        image=np.moveaxis(catchment_on_terrain, 0, -1),
        bounds=bbox,
        opacity=0.7,
        interactive=True,
        cross_origin=False,
        zindex=1,
    )
    # folium.Popup("Message").add_to(img)
    # Imposed Stream on the Terrain

    img = folium.raster_layers.ImageOverlay(
        name="Stream",
        image=np.moveaxis(array, 0, -1),
        bounds=bbox,
        opacity=0.3,
        interactive=True,
        cross_origin=False,
        zindex=1,
    )
    img.add_to(n)
    terrain_img.add_to(n)
    folium.LayerControl().add_to(n)
    folium_static(n)
else: pass

