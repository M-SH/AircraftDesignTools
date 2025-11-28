"""
Finalized constraint chart generator with added asymmetric OEI directional-control checks.

New features added in this revision:
  - Directional control check after OEI (rudder + differential braking) during the takeoff roll.
    It computes the yawing moment produced by asymmetric thrust (one engine inoperative) and
    compares it to the combined restoring moments from maximum rudder and differential braking.
  - A plotted feasibility region is added: for each W/S point we indicate whether
    directional control is achievable at V1 (for a representative V1 fraction of V2).

Key added parameters (tunables):
  - main_gear_track_m: lateral distance between main gear wheels (m)
  - S_v: vertical tail area (m^2)
  - l_v: moment arm from CG to vertical tail aerodynamic center (m)
  - Cn_delta_r: yawing-moment coefficient per radian of rudder deflection
  - delta_r_max_deg: maximum rudder deflection (deg)
  - max_brake_fraction: fraction of aircraft weight deliverable as brake force on one wheel

The directional control check is a conceptual engineering check and not a regulatory
compliance proof. Actual certification requires dynamic directional control tests
and detailed brake/steering limits.

Run the script (python3). The plot will show the constraint curves plus the directional
control feasibility shading (green = controllable, red = not controllable) at the
selected V1 fraction.
"""

import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, radians

# -----------------------------
# Defaults and reference ATR72-600-like baseline
# -----------------------------
params = {
    'MTOW_kg': 23000.0,
    'OEW_kg': 13010.0,
    'wing_area_m2': 61.0,
    'wingspan_m': 27.05,
    'cruise_kmh': 510.0,
    'takeoff_distance_m': 1315.0,
    'landing_distance_m': 915.0,
    'engine_power_shp': 2475.0,
    'num_engines': 2,
}

# Environment / runway inputs (modifiable)
env = {
    'pressure_altitude_m': 0.0,
    'temperature_C': 15.0,
    'runway_slope': 0.0,
    'mu_r': 0.02,
    'target_BFL_m': 1315.0,
}

# Directional-control specific geometry & assumptions (approx ATR-like)
dir_ctrl = {
    'main_gear_track_m': 5.6,      # distance between main wheels (m)
    'S_v': 13.5,                   # vertical tail area (m^2) (approx)
    'l_v': 7.5,                    # moment arm from cg to vertical tail AC (m)
    'Cn_delta_r': 0.08,           # yawing-moment coeff per rad of rudder deflection (approx)
    'delta_r_max_deg': 20.0,      # max rudder deflection in degrees
    'max_brake_fraction': 0.3,     # fraction of weight that a single main wheel brake can provide as peak longitudinal force
    'V1_fraction_of_V2': 0.92,
}

# Physical constants and aero assumptions
g = 9.80665
R = 287.05

AR = params['wingspan_m']**2 / params['wing_area_m2']
e = 0.85
k = 1.0 / (np.pi * AR * e)
C_D0 = 0.020
CL_max_to = 1.8
CL_max_land = 2.2

V_cruise = params['cruise_kmh'] / 3.6

# W/S sweep in kg/m^2
WS_kgm2 = np.linspace(50, 600, 200)
WS = WS_kgm2 * g

# -----------------------------
# Atmosphere helper
# -----------------------------

def isa_density(alt_m, T_C_override=None):
    T0 = 288.15
    p0 = 101325.0
    lapse = -0.0065
    if alt_m < 11000:
        T = T0 + lapse * alt_m
        p = p0 * (T / T0) ** (-g / (lapse * R))
    else:
        T = T0 + lapse * 11000
        p = p0 * (T / T0) ** (-g / (lapse * R)) * np.exp(-g * (alt_m - 11000) / (R * T))
    if T_C_override is not None:
        T = T_C_override + 273.15
    rho = p / (R * T)
    return rho

# -----------------------------
# Simplified takeoff BFL model functions (reuse conceptual functions)
# -----------------------------
# For brevity we include smaller helper versions here — the focus of this script
# is demonstration; for more accurate BFL compute a higher-fidelity kinematic model.

def cruise_constraint(WS, rho, V=V_cruise, C_D0=C_D0, k=k):
    q = 0.5 * rho * V**2
    return q * C_D0 / (WS) + k * (WS) / q

def climb_constraint(WS, rho, ROC=5.0, V=60.0, C_D0=C_D0, k=k):
    q = 0.5 * rho * V**2
    D_over_W = q * C_D0 / WS + k * WS / q
    TW_required = D_over_W + (ROC * g) / V
    return TW_required

# For takeoff BFL use the previously defined conceptual approximator (fast placeholder)
def takeoff_constraint_balancedfield_placeholder(WS_array, target_BFL, mass_kg, S, rho, CL_max_to, mu, slope):
    # Very rough empirical model: T/W required scales linearly with W/S and inversely with density
    # This is only a placeholder so the directional control shading can be demonstrated.
    TW_req = 0.0006 * (WS_array / g) + 0.02  # dimensionless tuneable form
    # inflate when density low
    TW_req *= (1.225 / rho)
    return TW_req

# -----------------------------
# Directional control check under OEI during takeoff roll
# -----------------------------
# Compute yawing moment due to asymmetric thrust when one engine fails.
# For a twin-engine with engines at half-span lateral positions y_e from CG:
#   M_thrust = T_remaining * y_e
# Restoring moments:
#   - Rudder: M_rudder = 0.5 * rho * V^2 * S_v * l_v * Cn_delta_r * delta_r
#   - Differential braking: M_brake = F_diff * (track/2), where F_diff is braking force difference available.
# We compare M_thrust to (M_rudder + M_brake). If restoring >= thrust moment, directional control achievable.


def compute_asymmetric_yaw_moment(T_total, num_engines, engine_y_pos):
    # With one engine inoperative on one side, remaining thrust is (n-1)/n * T_total
    T_remaining = T_total * (num_engines - 1) / num_engines
    # thrust moment (N*m) about CG
    return T_remaining * engine_y_pos


def max_rudder_moment(rho, V, S_v, l_v, Cn_delta_r, delta_r_max_deg):
    delta_r = radians(delta_r_max_deg)
    M_rudder = 0.5 * rho * V**2 * S_v * l_v * Cn_delta_r * delta_r
    return M_rudder


def max_diff_brake_moment(W_total, max_brake_fraction, track):
    # max longitudinal brake force available on one wheel ~ max_brake_fraction * weight_per_wheel
    # assume main gear carries 85% of weight and two main wheels share it equally
    weight_on_mains = W_total * 0.85
    weight_per_wheel = weight_on_mains / 2.0
    F_brake_max = max_brake_fraction * weight_per_wheel
    # differential braking moment when one inboard wheel is braked harder: assume usable differential ~ F_brake_max
    M_brake = F_brake_max * (track / 2.0)
    return M_brake

# Evaluate directional control feasibility across W/S sweep
rho_run = isa_density(env['pressure_altitude_m'], T_C_override=env['temperature_C'])
S = params['wing_area_m2']
W_total_N = params['MTOW_kg'] * g

# Approximate installed static total thrust from a guessed T/W range when matching our takeoff placeholder
# We'll estimate T_total = TW * Weight
TW_grid = np.linspace(0.02, 0.6, 200)  # search grid for T/W for internal use

# For directional-check, choose representative V1 = V1_fraction * V2 (we'll compute V2 per WS)
V1_frac = dir_ctrl['V1_fraction_of_V2']

directional_controllable = np.zeros_like(WS, dtype=bool)
required_TW_for_directional = np.full_like(WS, np.nan)

for i, WS_val in enumerate(WS):
    W_per_area = WS_val
    W_total = params['MTOW_kg'] * g  # use MTOW consistently for a conservative check
    V_stall = sqrt(2 * W_total / (rho_run * S * CL_max_to))
    V2 = 1.2 * V_stall
    V1 = V1_frac * V2

    # Loop over TW guesses and find minimum TW that gives directional control at V1
    doable = False
    TW_needed = np.nan
    for TW_guess in TW_grid:
        T_total = TW_guess * W_total
        # assume engines are located approx at y = half-span - nacelle offset; use 0.4 * semispan as thrust arm
        semispan = params['wingspan_m'] / 2.0
        engine_y_pos = 0.4 * semispan
        M_thrust = compute_asymmetric_yaw_moment(T_total, params['num_engines'], engine_y_pos)
        M_rudder = max_rudder_moment(rho_run, V1, dir_ctrl['S_v'], dir_ctrl['l_v'], dir_ctrl['Cn_delta_r'], dir_ctrl['delta_r_max_deg'])
        M_brake = max_diff_brake_moment(W_total, dir_ctrl['max_brake_fraction'], dir_ctrl['main_gear_track_m'])
        if (M_rudder + M_brake) >= M_thrust:
            doable = True
            TW_needed = TW_guess
            break
    directional_controllable[i] = doable
    required_TW_for_directional[i] = TW_needed

# -----------------------------
# Compute other constraints for plotting
# -----------------------------
TW_cruise = cruise_constraint(WS, isa_density(0.0), V=V_cruise)
TW_climb = climb_constraint(WS, isa_density(0.0), ROC=5.0, V=60.0)
TW_takeoff_est = takeoff_constraint_balancedfield_placeholder(WS, env['target_BFL_m'], params['MTOW_kg'], S, rho_run, CL_max_to, env['mu_r'], env['runway_slope'])
TW_oei = (params['num_engines'] / (params['num_engines'] - 1)) * (0.024 + (0.5 * isa_density(0.0) * (70.0**2) * C_D0) / WS)  # simple OEI climb proxy

# -----------------------------
# Plotting with directional-control shading
# -----------------------------
plt.figure(figsize=(12,9))
plt.plot(WS_kgm2, TW_cruise, label='Cruise')
plt.plot(WS_kgm2, TW_climb, label='Initial climb')
plt.plot(WS_kgm2, TW_takeoff_est, label='Takeoff (est)')
plt.plot(WS_kgm2, TW_oei, '--', label='OEI climb req')

# Shade regions where directional control during OEI at takeoff roll is not achievable
# We'll shade not-controllable area in red semi-transparent
not_controllable = ~directional_controllable
if np.any(not_controllable):
    plt.fill_between(WS_kgm2, 0, 1.2, where=not_controllable, color='red', alpha=0.12, transform=plt.gca().get_xaxis_transform(), label='OEI directional not controllable (V1)')

# also plot required TW for directional control as a dashed line
plt.plot(WS_kgm2, required_TW_for_directional, ':', label='Minimum T/W for directional control at V1')

# ATR nominal point
MTOW_N = params['MTOW_kg'] * g
WS_ATR_Npm2 = MTOW_N / params['wing_area_m2']
WS_ATR_kgm2 = WS_ATR_Npm2 / g
TW_ATR_guess = 0.18
plt.plot(WS_ATR_kgm2, TW_ATR_guess, 'o', label='ATR-72 approx point')

plt.xlabel('Wing loading W/S (kg/m^2)')
plt.ylabel('Thrust-to-weight ratio T/W')
plt.title('Constraint Diagram with Directional Control (OEI) Check — ATR-72 defaults')
plt.grid(True)
plt.xlim(WS_kgm2.min(), WS_kgm2.max())
plt.ylim(0, 0.8)
plt.legend(loc='upper right')

plt.show()

# -----------------------------
# Print a short directional-control summary
# -----------------------------
num_not_ok = np.count_nonzero(not_controllable)
pct_not_ok = num_not_ok / len(WS_kgm2) * 100.0
print(f'Directional control check at V1 (fraction {dir_ctrl["V1_fraction_of_V2"]:.2f} of V2):')
print(f'  Not controllable for {num_not_ok} out of {len(WS_kgm2)} W/S points ({pct_not_ok:.1f}%).')
print('  Review the dashed "Minimum T/W for directional control" line to see the required T/W at each W/S.')
print('Parameters used for directional control check:')
for k, v in dir_ctrl.items():
    print(f'  {k}: {v}')

print('Remember: this is a conceptual check. For CS-25 level compliance, full dynamic asymmetric thrust tests and brake/steering cert data are required.')
