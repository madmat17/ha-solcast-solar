from __future__ import annotations
import math
from dataclasses import dataclass
from typing import List, Tuple, Dict

# ----- Domain model (everything lives in the roof plane) ----------------------

@dataclass
class Dormer:
    # Dimensions in meters (footprint on roof plane)
    width: float           # along ridge (x-axis, left(-) to right(+))
    depth: float           # downslope (y-axis, ridge->eave)
    height: float          # dormer height above roof plane (meters)
    x_center: float        # dormer center x in roof coordinates (m, 0 at dormer midline)
    y_front: float         # y position of dormer front (top/up-slope edge), meters


@dataclass
class PanelArea:
    # Axis-aligned rectangles on roof plane (relative to dormer origin/axes)
    # x spans along ridge; y spans down-slope (same axes as above).
    x0: float; y0: float
    x1: float; y1: float

    @property
    def area(self) -> float:
        return max(0.0, self.x1 - self.x0) * max(0.0, self.y1 - self.y0)


@dataclass
class Roof:
    azimuth_deg: float  # 0=N, 90=E, 180=S, 270=W
    tilt_deg: float     # 0 = flat, 90 = vertical


@dataclass
class SiteGeometry:
    roof: Roof
    dormer: Dormer
    panels_left: PanelArea
    panels_right: PanelArea


# ----- Vector helpers ---------------------------------------------------------

def _deg2rad(d: float) -> float:
    return d * math.pi / 180.0

def _sun_vector(az_deg: float, el_deg: float) -> Tuple[float, float, float]:
    """Unit vector pointing from the surface toward the sun in world coords.
    azimuth: 0=N, 90=E; elevation above horizon."""
    az = _deg2rad(az_deg)
    el = _deg2rad(el_deg)
    x = math.cos(el) * math.sin(az)     # East(+)
    y = math.cos(el) * math.cos(az)     # North(+)
    z = math.sin(el)                    # Up(+)
    return (x, y, z)

def _roof_basis(az_deg: float, tilt_deg: float) -> Tuple[Tuple[float,float,float], Tuple[float,float,float], Tuple[float,float,float]]:
    """Return orthonormal basis for roof plane:
       ex: along ridge (left/right), ey: down-slope (ridge->eave), n: roof normal (pointing out of plane)."""
    # Roof "down-slope" direction is the projection of the roof aspect (azimuth) into the horizontal plane.
    # Roof faces azimuth 'az_deg' (its normal projected azimuth). Let 'aspect' be toward maximum downslope.
    az = _deg2rad(az_deg)
    # ey_horiz points 'downslope' in horizontal plane: from ridge toward eave
    ey_h = (math.sin(az), math.cos(az), 0.0)  # East, North, Up
    # roof normal tilted by 'tilt' from vertical toward -ey_h
    tilt = _deg2rad(tilt_deg)
    n = (
        -math.sin(tilt) * ey_h[0],
        -math.sin(tilt) * ey_h[1],
        math.cos(tilt)
    )
    # ex along ridge = n x ey (right-hand rule)
    # first get ey as unit vector on plane (combine horiz and small up)
    # The true ey (on plane) is perpendicular to ridge and on plane:
    # Choose ey to be "down-slope" unit on plane
    # Start with horizontal ey_h and rotate down by tilt around ex; simpler: make ey perpendicular to n and to ex.
    # We'll pick temporary ex' as (-ey_h_y, ey_h_x, 0) (horiz ridge), then orthonormalize.
    ex_h = (-ey_h[1], ey_h[0], 0.0)
    # Orthonormalize ex against n
    def _norm(v):
        l = math.sqrt(sum(c*c for c in v))
        return (v[0]/l, v[1]/l, v[2]/l)
    def _dot(a,b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
    def _sub(a,b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    ex = _norm(_sub(ex_h, tuple(_dot(ex_h,n)*c for c in n)))
    # Now ey = n x ex (downslope)
    ey = _norm((n[1]*ex[2]-n[2]*ex[1], n[2]*ex[0]-n[0]*ex[2], n[0]*ex[1]-n[1]*ex[0]))
    return ex, ey, _norm(n)

def _project_to_roof_plane(v: Tuple[float,float,float], ex, ey, n) -> Tuple[float,float,float,float,float]:
    """Return: dot with n, projected vector components (vx, vy) on plane, and its magnitude."""
    # v dot n: cosine between sun vector and roof normal
    vdotn = v[0]*n[0] + v[1]*n[1] + v[2]*n[2]
    # Projection onto plane: v_proj = v - (v·n) n
    vx3 = v[0] - vdotn*n[0]
    vy3 = v[1] - vdotn*n[1]
    vz3 = v[2] - vdotn*n[2]
    # Components in roof basis
    vx = vx3*ex[0]+vy3*ex[1]+vz3*ex[2]  # along ridge
    vy = vx3*ey[0]+vy3*ey[1]+vz3*ey[2]  # downslope
    mag = math.hypot(vx, vy)
    return vdotn, vx, vy, mag, vdotn

# ----- Shadow footprint & overlap ---------------------------------------------

def _rect_overlap(a: PanelArea, b: PanelArea) -> float:
    x_overlap = max(0.0, min(a.x1, b.x1) - max(a.x0, b.x0))
    y_overlap = max(0.0, min(a.y1, b.y1) - max(a.y0, b.y0))
    return x_overlap * y_overlap

def _shadow_rect(geom: SiteGeometry, sun_az: float, sun_el: float) -> PanelArea | None:
    """Approximate the dormer shadow on roof plane as an axis-aligned rectangle in roof coords.
       Width is dormer width. Length is determined by sun vector projected onto the plane.
       If sun behind roof or below horizon, return None."""
    if sun_el <= 0.0:
        return None

    ex, ey, n = _roof_basis(geom.roof.azimuth_deg, geom.roof.tilt_deg)
    sv = _sun_vector(sun_az, sun_el)
    # If sun is 'behind' the roof, v·n <= 0 -> no direct sun; shadow not applicable for generation (roof self-shaded)
    vdotn, vx, vy, vmag, _ = _project_to_roof_plane(sv, ex, ey, n)
    if vdotn <= 0.0 or vmag == 0.0:
        return None

    # Shadow length along -projection direction. For an object offset 'height' along +n,
    # length L on plane = height * |v_proj| / (v · n)
    L = geom.dormer.height * (vmag / vdotn)

    # Shadow is cast in direction opposite to projected sun direction
    dx = -vx / vmag if vmag > 0 else 0.0
    dy = -vy / vmag if vmag > 0 else 0.0

    # Dormer footprint front/top edge at (x_center, y_front); extend rectangle of width and depth
    half_w = geom.dormer.width / 2.0
    dormer_rect = PanelArea(
        x0=geom.dormer.x_center - half_w,
        y0=geom.dormer.y_front,
        x1=geom.dormer.x_center + half_w,
        y1=geom.dormer.y_front + geom.dormer.depth,
    )

    # Build a shadow rectangle starting at the front edge of the dormer and extending L
    # Use axis-aligned bounding box of the swept area for simplicity.
    # Front line mid-point:
    x_mid = geom.dormer.x_center
    y_front = geom.dormer.y_front
    # Shadow bbox
    x0 = min(dormer_rect.x0, dormer_rect.x1, x_mid + dx*L - half_w, x_mid + dx*L + half_w)
    x1 = max(dormer_rect.x0, dormer_rect.x1, x_mid + dx*L - half_w, x_mid + dx*L + half_w)
    y0 = min(y_front, y_front + dy*L)
    y1 = max(y_front, y_front + dy*L)

    return PanelArea(x0=x0, y0=y0, x1=x1, y1=y1)

def compute_weighted_dampening(
    geom: SiteGeometry,
    sun_positions: List[Tuple[float, float]],  # list[(azimuth_deg, elevation_deg)] for each interval
) -> List[float]:
    """Return [0..1] dampening factors per interval (1 = no dampening).
       Weighted by left/right panel areas."""
    total_area = geom.panels_left.area + geom.panels_right.area
    if total_area <= 0.0:
        return [1.0 for _ in sun_positions]

    factors: List[float] = []
    for az, el in sun_positions:
        sh = _shadow_rect(geom, az, el)
        if sh is None:
            factors.append(1.0)
            continue

        # Overlaps
        ov_l = _rect_overlap(sh, geom.panels_left)
        ov_r = _rect_overlap(sh, geom.panels_right)

        # Clip to panel areas (numeric noise guard)
        ov_l = min(ov_l, geom.panels_left.area)
        ov_r = min(ov_r, geom.panels_right.area)

        # Per-side uncovered fraction
        frac_l = 1.0 - (ov_l / geom.panels_left.area if geom.panels_left.area > 0 else 0.0)
        frac_r = 1.0 - (ov_r / geom.panels_right.area if geom.panels_right.area > 0 else 0.0)

        # Weighted by area
        weighted = ((frac_l * geom.panels_left.area) + (frac_r * geom.panels_right.area)) / total_area
        # Damping factor is the proportion of generation not eclipsed
        # clamp into [0,1]
        factors.append(max(0.0, min(1.0, weighted)))
    return factors
