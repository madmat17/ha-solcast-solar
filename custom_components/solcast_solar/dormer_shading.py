import math

def sun_position(latitude, longitude, datetime_obj):
    """Returns sun azimuth and elevation at the given location and time."""
    from astral import LocationInfo
    from astral.sun import sun
    location = LocationInfo(latitude=latitude, longitude=longitude)
    s = sun(location.observer, date=datetime_obj)
    return s['azimuth'], s['elevation']

def project_shadow_on_roof(dormer, roof_tilt, sun_azimuth, sun_elevation, roof_azimuth):
    """
    Projects the dormer shadow onto the roof. Returns shadow start/end coordinates (meters from roof left edge).
    dormer: dict with 'width', 'height', 'depth', 'distance_from_left_edge', 'distance_from_eave'
    roof_tilt: degrees
    sun_azimuth: degrees
    sun_elevation: degrees
    roof_azimuth: degrees
    """
    # Calculate shadow length on tilted roof
    tilt_rad = math.radians(roof_tilt)
    sun_el_rad = math.radians(sun_elevation)
    # Effective sun elevation over roof plane
    effective_sun_el = math.asin(math.sin(sun_el_rad) * math.cos(tilt_rad) +
                                 math.cos(sun_el_rad) * math.sin(tilt_rad) * math.cos(math.radians(sun_azimuth - roof_azimuth)))
    effective_sun_el = max(effective_sun_el, math.radians(5))  # Avoid very small angles
    shadow_length = dormer['height'] / math.tan(effective_sun_el)

    # Shadow direction: project along sun azimuth minus roof azimuth
    direction = math.radians(sun_azimuth - roof_azimuth)
    # Shadow projected onto roof width
    lateral_offset = shadow_length * math.sin(direction)
    # Shadow projected along roof from dormer front edge
    shadow_start = dormer['distance_from_left_edge']
    shadow_end = shadow_start + dormer['width'] + lateral_offset

    # Clamp to roof boundaries if needed
    return min(shadow_start, shadow_end), max(shadow_start, shadow_end)

def compute_pv_shadow_fraction(shadow_start, shadow_end, pv_left, pv_right):
    """
    Returns fraction of shadowed area on each PV field.
    pv_left/pv_right: dict with 'width', 'height', 'distance_from_left_edge'
    Assumes panels are rectangles aligned along roof width.
    """
    # PV left
    pv_left_start = pv_left['distance_from_left_edge']
    pv_left_end = pv_left_start + pv_left['width']
    overlap_left = max(0, min(shadow_end, pv_left_end) - max(shadow_start, pv_left_start))
    frac_left = overlap_left / pv_left['width'] if pv_left['width'] > 0 else 0

    # PV right
    pv_right_start = pv_right['distance_from_left_edge']
    pv_right_end = pv_right_start + pv_right['width']
    overlap_right = max(0, min(shadow_end, pv_right_end) - max(shadow_start, pv_right_start))
    frac_right = overlap_right / pv_right['width'] if pv_right['width'] > 0 else 0

    return frac_left, frac_right

def calculate_dynamic_damping(dormer, pv_left, pv_right, roof_azimuth, roof_tilt, latitude, longitude, times):
    """
    dormer: dict with 'width', 'height', 'depth', 'distance_from_left_edge', 'distance_from_eave'
    pv_left: dict with 'width', 'height', 'distance_from_left_edge'
    pv_right: dict with 'width', 'height', 'distance_from_left_edge'
    roof_azimuth: degrees
    roof_tilt: degrees
    latitude, longitude: site coordinates
    times: list of datetime objects
    """
    damping_factors = []
    pv_left_area = pv_left['width'] * pv_left['height']
    pv_right_area = pv_right['width'] * pv_right['height']
    total_area = pv_left_area + pv_right_area
    for dt in times:
        sun_az, sun_el = sun_position(latitude, longitude, dt)
        shadow_start, shadow_end = project_shadow_on_roof(
            dormer, roof_tilt, sun_az, sun_el, roof_azimuth
        )
        frac_left, frac_right = compute_pv_shadow_fraction(
            shadow_start, shadow_end, pv_left, pv_right
        )
        weighted_damping = ((frac_left * pv_left_area) + (frac_right * pv_right_area)) / total_area
        damping_factors.append(weighted_damping)
    avg_damping = sum(damping_factors) / len(damping_factors)
    return avg_damping
