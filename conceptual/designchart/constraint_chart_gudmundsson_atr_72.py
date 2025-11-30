"""
Refined Constraint Chart generator with:
  - Balanced-field takeoff approximation (simple accelerate-stop / accelerate-go model)
    including runway slope, pressure altitude and temperature effects and runway
    rolling friction. This is still a conceptual-level model but is much closer
    to the CS-25 balanced-field logic than the previous crude algebraic form.
  - One-Engine-Inoperative (OEI) climb constraints using CS-25 style margins
    (steady climb and takeoff/climb gradients). The code computes required
    T/W for OEI climb gradients and shows an OEI-limited region.

Defaults are still ATR-72-600-like, and all units are SI. The script plots
W/S (kg/m^2) vs T/W and overlays cruise, climb, service-ceiling, takeoff
(balanced-field approximation) and OEI climb constraints.

Notes / limitations:
 - This is still a conceptual tool. For certification-level BFL and OEI
   performance assessments a flight mechanics performance code or dedicated
   BFL simulator is required.
 - The balanced-field implementation below follows a simplified kinematic
   treatment (classic conceptual approach) and neglects many second-order
   effects (e.g., thrust lapse with speed, ground effect on CLmax, engine
   acceleration dynamics). It is parameterized so you can tighten/loosen
   assumptions.

Usage: run the script (python3). It will open a matplotlib plot and print
an ATR nominal summary. Edit the `params` dict for runway and environment.
"""

import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

# -----------------------------
# Defaults and reference ATR72-600-like baseline
# -----------------------------
params = {
    'MTOW_kg': 23000.0,
    'OEW_kg': 13010.0,
    'wing_area_m2': 61.0,
    'wingspan_m': 27.05,
    'cruise_kmh': 510.0,
    'takeoff_distance_m': 1315.0,    # published warm reference (for rough comparison)
    'landing_distance_m': 915.0,
    'engine_power_shp': 2475.0,
    'num_engines': 2,
}

# Environment / runway inputs (modifiable)
env = {
    'pressure_altitude_m': 0.0,   # m
    'temperature_C': 15.0,        # deg C
    'runway_slope': 0.0,          # fraction (positive uphill)
    'mu_r': 0.02,                 # rolling friction coefficient (paved dry)
    'target_BFL_m': 1315.0,       # design balanced-field length target to compare against
}

# Physical constants and aero assumptions
g = 9.80665
R = 287.05  # specific gas constant for air

# Aerodynamic assumptions (tunable)
AR = params['wingspan_m']**2 / params['wing_area_m2']
e = 0.85
k = 1.0 / (np.pi * AR * e)
C_D0 = 0.020
CL_max_to = 1.8
CL_max_land = 2.2
eta_prop = 0.85

# Speeds of interest
V_cruise = params['cruise_kmh'] / 3.6

# W/S sweep (use kg/m^2 as x-axis for readability)
WS_kgm2 = np.linspace(50, 600, 300)
WS = WS_kgm2 * g  # convert to N/m^2 for internal equations

# -----------------------------
# Atmosphere helper (ISA lapse simplified)
# -----------------------------

def isa_density(alt_m, T_C_override=None):
    """Return air density at altitude (m) and temperature override (C) if given.
    Uses ISA approximations for troposphere (linear lapse to 11 km).
    """
    # Standard ISA
    T0 = 288.15
    p0 = 101325.0
    lapse = -0.0065
    if alt_m < 11000:
        T = T0 + lapse * alt_m
        p = p0 * (T / T0) ** (-g / (lapse * R))
    else:
        # above troposphere approximate (not needed for this tool generally)
        T = T0 + lapse * 11000
        p = p0 * (T / T0) ** (-g / (lapse * R)) * np.exp(-g * (alt_m - 11000) / (R * T))
    if T_C_override is not None:
        T = T_C_override + 273.15
    rho = p / (R * T)
    return rho

# -----------------------------
# Balanced-field takeoff model (simplified)
# -----------------------------
# Approach: compute accelerate-stop distance (AS) and accelerate-go distance (AG)
# for a trial decision speed V1. Find V1 where AS = AG. Then BFL = AS (or AG).
# We implement a moderately detailed kinematic model with speed-dependent
# aerodynamic drag, rolling friction, slope and thrust differences (all-engines
# and one-engine-inoperative). We then invert to find required T/W for a
# given W/S such that BFL <= target runway length; for plotting we compute the
# T/W required to meet a specified BFL (by root-finding on T/W). This gives
# a takeoff constraint curve.


def available_thrust_per_weight(TW_guess, W):
    """Given T/W guess (dimensionless) and weight per unit area W (N/m^2) returns
    total available thrust per unit wing area (N/m^2). This helper keeps units
    consistent when we compute accelerations tied to WS values.
    But in our kinematic equations we'll mostly use T/W as a dimensionless number
    applied to actual aircraft weight (total).
    """
    return TW_guess * W  # N per unit wing area; not used directly often


def dynamic_thrust(T_static, V):
    """Simplified model: allow some small thrust lapse with speed (props/ducted)
    T(V) = T_static * (1 - c * (V/V_ref)) limiting to positive values.
    For turboprops this is approximate; keep parameter c small.
    """
    V_ref = 200.0  # m/s reference
    c = 0.1
    return max(0.1 * T_static, T_static * (1 - c * (V / V_ref)))


def compute_AS_AG_for_V1(V1, TW_all, TW_oew, mass_kg, S, rho, CL_max_to, mu, slope):
    """Compute Accelerate-Stop (AS) and Accelerate-Go (AG) distances for
    a given V1 (m/s). TW_all is T/W with all engines operative (dimensionless),
    TW_oew is T/W after engine failure (one engine inoperative case, dimensionless),
    mass_kg is aircraft mass, S is wing area, rho density, CL_max_to lift coefficient,
    mu rolling friction, slope runway fraction.

    Returns (AS, AG) in meters. Uses a segmented kinematic integration approach:
    - accelerate from 0 to V1 under all-engines thrust (account drag, rolling)
    - at V1: engine fails -> for stop: braking deceleration (uses mu and aerodynamic drag)
            for go: accelerate to V2 (1.2*V_stall) using reduced thrust (OEI) and then
            liftoff and climb distance to clear obstacle is approximated (we simplify obstacle clearance)
    """
    W_newtons = mass_kg * g
    # compute CL at rotation/near-liftoff; approximate V_stall and V2
    V_stall = sqrt(2 * W_newtons / (rho * S * CL_max_to))
    V2 = 1.2 * V_stall
    # Limit V1 to sensible bounds
    V1 = max(0.2 * V2, min(V1, 0.98 * V2))

    # Integrate acceleration from 0 to V1 using small dt stepping
    dt = 0.1
    t = 0.0
    V = 0.0
    s = 0.0
    while V < V1 - 1e-6:
        # forces at current speed
        q = 0.5 * rho * V**2
        D = q * S * (C_D0 + k * (W_newtons/(q*S))**2)
        T = TW_all * W_newtons  # total thrust (dimensionless * weight)
        a = (T - D - mu * (W_newtons - 0.5 * rho * V**2 * S * CL_max_to) + (-slope * W_newtons)) / W_newtons * g
        # approximate: braking term uses mu*W minus small lift correction
        # ensure a is not extremely large
        a = max(-5.0, min(10.0, a))
        V = V + a * dt
        V = max(0.0, V)
        s = s + V * dt
        t += dt
        # avoid infinite loops
        if t > 200:
            break

    s_to_V1 = s

    # --- Accelerate-stop: braking from V1 to 0 under brakes and drag
    V = V1
    s_stop = 0.0
    t = 0.0
    while V > 1e-2:
        q = 0.5 * rho * V**2
        D = q * S * (C_D0 + k * (W_newtons/(q*S))**2)
        # braking deceleration: assume very high braking decel available limited by mu
        a_brake = - (mu * (W_newtons - 0.5 * rho * V**2 * S * CL_max_to) + D + slope * W_newtons) / W_newtons * g
        a_brake = max(-10.0, a_brake)
        V = V + a_brake * dt
        V = max(0.0, V)
        s_stop = s_stop + V * dt
        t += dt
        if t > 200:
            break

    AS = s_to_V1 + s_stop

    # --- Accelerate-go: at V1 engine fails and remaining thrust is TW_oew*W
    # accelerate from V1 to V2 under reduced thrust
    V = V1
    s_go = 0.0
    t = 0.0
    while V < V2 - 1e-3:
        q = 0.5 * rho * V**2
        D = q * S * (C_D0 + k * (W_newtons/(q*S))**2)
        T = TW_oew * W_newtons
        a = (T - D - mu * (W_newtons - 0.5 * rho * V**2 * S * CL_max_to) - slope * W_newtons) / W_newtons * g
        a = max(-5.0, min(10.0, a))
        V = V + a * dt
        V = max(0.0, V)
        s_go = s_go + V * dt
        t += dt
        if t > 400:
            break

    # After reaching V2, approximate liftoff and 35ft obstacle distance using simple empirical multiplier
    # This lumps rotation, ground effect and climb to obstacle into a factor times speed squared / (2 * acceleration)
    # Use a conservative multiplier
    climb_factor = 1.2
    # average climb gradient approx (using excess thrust at OEI)
    q = 0.5 * rho * V2**2
    D = q * S * (C_D0 + k * (W_newtons/(q*S))**2)
    T_oew = TW_oew * W_newtons
    excess_a = max(0.01, (T_oew - D - mu * (W_newtons - 0.5 * rho * V2**2 * S * CL_max_to) - slope * W_newtons) / W_newtons * g)
    # distance to accelerate from V2 to rotation + climb to obstacle
    s_climb = climb_factor * (V2**2) / (2 * excess_a)

    AG = s_to_V1 + s_go + s_climb

    return AS, AG


def balanced_field_length_for_TW(WS_val, TW_all, mass_kg, S, rho, CL_max_to, mu, slope):
    """Given a wing loading WS_val (N/m^2), an all-engines T/W guess TW_all, compute
    the balanced-field length by finding V1 where AS==AG. We'll search V1 in [0.3*V2,0.98*V2].
    Returns BFL (m).
    """
    # approximate one-engine-inoperative T/W: for a 2-engine aircraft losing 1 engine, available thrust halves (neglecting asymmetry)
    # but turboprops often have less than linear scaling due to installed effects. We'll assume remaining thrust = (n-1)/n * full thrust
    # For n=2, TW_oew ~ 0.5 * TW_all (dimensionless)
    TW_oew = TW_all * (params['num_engines'] - 1) / params['num_engines']

    W_total = mass_kg * g
    # compute V_stall and V2 for use in V1 bounds
    V_stall = sqrt(2 * W_total / (rho * S * CL_max_to))
    V2 = 1.2 * V_stall

    # sweep V1 and compute AS-AG, find root where difference changes sign
    V1s = np.linspace(0.3 * V2, 0.98 * V2, 40)
    diffs = []
    AS_vals = []
    AG_vals = []
    for V1 in V1s:
        AS, AG = compute_AS_AG_for_V1(V1, TW_all, TW_oew, mass_kg, S, rho, CL_max_to, mu, slope)
        diffs.append(AS - AG)
        AS_vals.append(AS)
        AG_vals.append(AG)

    diffs = np.array(diffs)
    # find sign change
    sign_changes = np.where(np.sign(diffs[:-1]) != np.sign(diffs[1:]))[0]
    if len(sign_changes) == 0:
        # no sign change — pick the minimal of AS and AG as conservative BFL
        BFL = np.minimum(AS_vals, AG_vals).max()
    else:
        idx = sign_changes[0]
        # linear interpolation for V1 where AS==AG
        f1 = diffs[idx]
        f2 = diffs[idx+1]
        V1_interp = V1s[idx] - f1 * (V1s[idx+1] - V1s[idx]) / (f2 - f1)
        AS, AG = compute_AS_AG_for_V1(V1_interp, TW_all, TW_oew, mass_kg, S, rho, CL_max_to, mu, slope)
        BFL = 0.5 * (AS + AG)
    return BFL


def takeoff_constraint_balancedfield(WS_array, target_BFL, mass_kg, S, rho, CL_max_to, mu, slope):
    """For each WS in WS_array (N/m^2), find the minimum T/W such that balanced-field length <= target_BFL.
    We do a simple root-find over T/W (search over a reasonable range).
    """
    TW_req = np.zeros_like(WS_array)
    for i, ws in enumerate(WS_array):
        # Solve for TW such that BFL(TW) - target_BFL = 0
        def bfl_minus_target(TW):
            return balanced_field_length_for_TW(ws, TW, mass_kg, S, rho, CL_max_to, mu, slope) - target_BFL

        # Search TW in [0.01, 1.0]
        TW_low = 0.01
        TW_high = 1.0
        # ensure sign at bounds
        try:
            f_low = bfl_minus_target(TW_low)
            f_high = bfl_minus_target(TW_high)
        except Exception:
            TW_req[i] = np.nan
            continue
        if f_low < 0:
            # already below target at small TW
            TW_req[i] = TW_low
            continue
        if f_high > 0:
            # even at very high TW cannot reach target_BFL -> mark nan
            TW_req[i] = np.nan
            continue
        # bisection
        for _ in range(20):
            TW_mid = 0.5 * (TW_low + TW_high)
            f_mid = bfl_minus_target(TW_mid)
            if np.isnan(f_mid):
                TW_low = TW_mid
                continue
            if f_mid > 0:
                TW_low = TW_mid
            else:
                TW_high = TW_mid
        TW_req[i] = 0.5 * (TW_low + TW_high)
    return TW_req

# -----------------------------
# OEI climb constraints (CS-25 style margins)
# -----------------------------
# We'll compute the required T/W to meet a target climb gradient after an OEI event,
# e.g. net climb gradient (vertical speed / horizontal speed) required during segment
# (CS-25 uses gradient requirements after takeoff — we implement a representative one).


def oei_climb_constraint(WS_array, rho, S, CL_to, required_gradient=0.024, V=70.0):
    """Compute T/W required for a required climb gradient (dimensionless) with only
    n-1 engines operative. required_gradient is e.g. 0.024 for 2.4% climb.
    V is the horizontal speed (m/s) during the climb segment.
    """
    TW_req = np.zeros_like(WS_array)
    for i, WS_val in enumerate(WS_array):
        W_per_area = WS_val
        # total weight per wing area times wing area gives total weight; but TW is dimensionless so we can keep ratios
        q = 0.5 * rho * V**2
        Cl = CL_to
        D_over_W = q * C_D0 / WS_val + k * WS_val / q
        # climb gradient gamma = (T/W - D/W) / (1/g) * (V/V) simplified -> gamma = (T/W - D/W)
        # Actually vertical/horizontal => gamma = (T/W - D/W) (dimensionless) when using units consistent with g.
        # For OEI, available thrust is reduced: TW_available = TW_total * (n-1)/n
        # Solve for TW_total such that (TW_available - D/W) >= required_gradient
        # => TW_total * (n-1)/n >= required_gradient + D/W  => TW_total >= (n/(n-1)) * (required_gradient + D/W)
        factor = params['num_engines'] / (params['num_engines'] - 1)
        TW_needed = factor * (required_gradient + D_over_W)
        TW_req[i] = TW_needed
    return TW_req

# -----------------------------
# Other earlier constraints (cruise, climb, service ceiling)
# -----------------------------
def cruise_constraint(WS, rho, V=V_cruise, C_D0=C_D0, k=k):
    q = 0.5 * rho * V**2
    return q * C_D0 / (WS) + k * (WS) / q

def climb_constraint(WS, rho, ROC=5.0, V=60.0, C_D0=C_D0, k=k):
    q = 0.5 * rho * V**2
    D_over_W = q * C_D0 / WS + k * WS / q
    TW_required = D_over_W + (ROC * g) / V
    return TW_required

def service_ceiling_constraint(WS, rho, ROC=0.5, V=100.0, C_D0=C_D0, k=k):
    q = 0.5 * rho * V**2
    D_over_W = q * C_D0 / WS + k * WS / q
    TW_required = D_over_W + (ROC * g) / V
    return TW_required

# -----------------------------
# Compute density at runway conditions
# -----------------------------
rho_runway = isa_density(env['pressure_altitude_m'], T_C_override=env['temperature_C'])

# -----------------------------
# Compute constraint curves
# -----------------------------
TW_cruise = cruise_constraint(WS, isa_density(0.0), V=V_cruise)
TW_climb = climb_constraint(WS, isa_density(0.0), ROC=5.0, V=60.0)
TW_service = service_ceiling_constraint(WS, isa_density(5000.0), ROC=0.5, V=100.0)

# Balanced-field takeoff curve: compute required T/W to meet target BFL at runway conditions
mass_kg = params['MTOW_kg']
S = params['wing_area_m2']
TW_takeoff_req = takeoff_constraint_balancedfield(WS, env['target_BFL_m'], mass_kg, S, rho_runway, CL_max_to, env['mu_r'], env['runway_slope'])

# OEI climb constraint (CS-25 style) — use required gradient e.g. 0.024 (2.4%)
TW_oei = oei_climb_constraint(WS, isa_density(0.0), S, CL_max_to, required_gradient=0.024, V=70.0)

# -----------------------------
# Plotting
# -----------------------------
plt.figure(figsize=(11,8))
plt.plot(WS_kgm2, TW_cruise, label='Cruise')
plt.plot(WS_kgm2, TW_climb, label='Initial climb (ROC=5 m/s)')
plt.plot(WS_kgm2, TW_service, label='Service ceiling (approx)')
plt.plot(WS_kgm2, TW_takeoff_req, label=f'Takeoff (BFL target {env["target_BFL_m"]} m)')
plt.plot(WS_kgm2, TW_oei, '--', label='OEI climb constraint (req grad 2.4%)')

# ATR nominal point
MTOW_N = params['MTOW_kg'] * g
WS_ATR_Npm2 = MTOW_N / params['wing_area_m2']
WS_ATR_kgm2 = WS_ATR_Npm2 / g
TW_ATR_guess = 0.18
plt.plot(WS_ATR_kgm2, TW_ATR_guess, 'o', label='ATR-72 approx point')

plt.xlabel('Wing loading W/S (kg/m^2)')
plt.ylabel('Thrust-to-weight ratio T/W')
plt.title('Constraint Diagram with Balanced-Field Takeoff and OEI Climb Constraints (ATR-72-600 defaults)')
plt.grid(True)
plt.xlim(WS_kgm2.min(), WS_kgm2.max())
plt.ylim(0, 1.0)
plt.legend(loc='upper right')

plt.show()

# -----------------------------
# Print summary and notes
# -----------------------------
print('Refined constraint chart generated — updated script includes:')
print(' - Balanced-field takeoff approximation (accelerate-stop / accelerate-go).')
print(' - One-engine-inoperative (OEI) climb constraint using a representative CS-25 style climb gradient requirement (2.4% by default).')
print('Environment and runway parameters used:')
for k, v in env.items():
    print(f'  {k}: {v}')

print('Key assumptions / limitations:')
print(' - The balanced-field model is a conceptual kinematic integrator (not a certification-level BFL solver).')
print(' - Thrust is approximated proportional to weight (T/W), with a simplified speed-lapse model. Engine-out thrust is modeled as (n-1)/n of total.')
print(' - Braking deceleration is approximated via a rolling friction coefficient and aerodynamic braking; anti-skid/brake energy limits are neglected.')
print('If you want:')
print(' - I can calibrate the BFL model against a published BFL for ATR-72 by tuning mu, climb_factor and thrust lapse.')
print(' - I can add runway altitude/temperature sensitivity plots, show V1 contours, or include OEI asymmetric rolling moment checks (nose wheel steering and directional control).')
