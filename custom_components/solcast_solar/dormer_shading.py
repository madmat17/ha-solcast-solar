import math
from astral import LocationInfo
from astral.sun import sun

def sun_position(latitude, longitude, datetime_obj):
    """Returns sun azimuth and elevation at the given location and time."""
    location = LocationInfo(latitude=latitude, longitude=longitude)
    s = sun(location.observer, date=datetime_obj)
    return s['azimuth'], s['elevation']

def project_shadow_on_roof(dormer, roof_tilt, sun_azimuth, sun_elevation, roof_azimuth):
    """
    Projects the dormer shadow onto the roof. Returns shadow start/end coordinates (meters from roof left edge).
    """
    # Effective sun elevation over roof plane
    tilt_rad = math.radians(roof_tilt)
    sun_el_rad = math.radians(sun_elevation)
    sun_az_rel = math.radians(sun_azimuth - roof_azimuth)
    effective_sun_el = math.asin(
        math.sin(sun_el_rad) * math.cos(tilt_rad) +
        math.cos(sun_el_rad) * math.sin(tilt_rad) * math.cos(sun_az_rel)
    )
    effective_sun_el = max(effective_sun_el, math.radians(5))
    shadow_length = dormer['height'] / math.tan(effective_sun_el)

    # Project shadow along roof's width axis
    direction = sun_az_rel
    lateral_offset = shadow_length * math.sin(direction)
    shadow_start = dormer['distance_from_left_edge']
    shadow_end = shadow_start + dormer['width'] + lateral_offset

    return min(shadow_start, shadow_end), max(shadow_start, shadow_end)

def compute_pv_shadow_fraction(shadow_start, shadow_end, pv_field):
    """
    Returns fraction of shadowed area on one PV field.
    """
    pv_start = pv_field['distance_from_left_edge']
    pv_end = pv_start + pv_field['width']
    overlap = max(0, min(shadow_end, pv_end) - max(shadow_start, pv_start))
    covered_area = overlap * pv_field['height']
    total_area = pv_field['width'] * pv_field['height']
    fraction = covered_area / total_area if total_area > 0 else 0
    return fraction

def calculate_dynamic_damping(dormer, pv_left, pv_right, roof_azimuth, roof_tilt, latitude, longitude, times):
    """
    Returns the average weighted damping factor over a day for the given geometry and sun positions.
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
        frac_left = compute_pv_shadow_fraction(shadow_start, shadow_end, pv_left)
        frac_right = compute_pv_shadow_fraction(shadow_start, shadow_end, pv_right)
        weighted_damping = ((frac_left * pv_left_area) + (frac_right * pv_right_area)) / total_area
        damping_factors.append(weighted_damping)
    avg_damping = sum(damping_factors) / len(damping_factors)
    return avg_damping
