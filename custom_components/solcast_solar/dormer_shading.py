import math
from astral.sun import sun

def calculate_dynamic_damping(dormer, pv_left, pv_right, roof_azimuth, roof_angle, times, location):
    """Calculate dynamic damping factors based on dormer geometry, PV panel measurements,
    roof orientation, and sun's position for the day."""
    damping_factors = []
    for dt in times:
        sun_pos = sun(location.observer, date=dt)
        sun_elev = sun_pos["elevation"]
        sun_az = sun_pos["azimuth"]

        # Geometry logic (simplified, refine as needed)
        shadow_length = dormer["height"] / math.tan(math.radians(max(sun_elev, 1)))
        shadow_offset = dormer["distance_from_left_edge"] + shadow_length * math.sin(math.radians(sun_az - roof_azimuth))
        pv_left_covered = max(0, min(pv_left["width"], shadow_offset)) / pv_left["width"]
        pv_right_covered = max(0, min(pv_right["width"], shadow_offset)) / pv_right["width"]
        total_area = pv_left["width"] * pv_left["height"] + pv_right["width"] * pv_right["height"]
        weighted_damping = (
            pv_left_covered * pv_left["width"] * pv_left["height"] +
            pv_right_covered * pv_right["width"] * pv_right["height"]
        ) / total_area
        damping_factors.append(weighted_damping)
    avg_damping = sum(damping_factors) / len(damping_factors)
    return avg_damping
