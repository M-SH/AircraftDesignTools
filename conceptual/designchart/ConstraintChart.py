# cs25_constraint_a320.py
# Run in any modern Python 3 environment with numpy and matplotlib installed.
import numpy as np
import matplotlib.pyplot as plt

# ----------------- User switches -----------------
use_imperial = False        # True -> show wing loading axis in lb/ft^2
save_png = False            # True -> save diagram to ./constraint_a320.png

# ----------------- Physical constants -----------------
g = 9.80665                # m/s^2
rho_sl = 1.225             # kg/m^3 sea-level ISA (approx)

# ----------------- Aerodynamic helpers -----------------
def q_dynamic(rho, V):
    return 0.5 * rho * V**2

def CL_for_WS_at_speed(ws, rho, V):
    # CL needed to support W/S at speed V: CL = 2*(W/S)/(rho*V^2)
    return 2.0 * ws / (rho * V**2)

def drag_over_weight(ws, rho, V, CD0, k):
    # D/W = (q * (CD0 + k*CL^2)) / (W/S)
    q = q_dynamic(rho, V)
    CL = CL_for_WS_at_speed(ws, rho, V)
    CD = CD0 + k * CL**2
    return (q * CD) / ws

def stall_speed(ws, rho, CL_max):
    # V_stall = sqrt( 2*(W/S) / (rho * CL_max) )
    return np.sqrt(2.0 * ws / (rho * CL_max))

# ----------------- Takeoff / BFL approximations -----------------
def takeoff_ground_roll(Vr, TW, mu_roll=0.02):
    # average acceleration a = g*(T/W - mu_roll), approximate ground roll: s = Vr^2/(2*a)
    a = g * (TW - mu_roll)
    # avoid negative/zero accelerations
    if a <= 0:
        return np.inf
    return Vr**2 / (2.0 * a)

def takeoff_airborne_distance(Vr, climb_angle_est=0.03, factor=1.2):
    # Simple energetic approx for airborne distance to clear obstacle:
    # s_air ~ factor * Vr^2 / (2*g) scaled by climb angle losses.
    # factor ~ 1.0..1.3 to capture drag/transition losses
    return factor * Vr**2 / (2.0 * g) * (1.0 + 1.0/climb_angle_est)

def stop_distance_from_speed(V, mu_brake=0.4, t_reaction=1.5):
    # s_reaction + braking distance
    s_reaction = V * t_reaction
    s_brake = V**2 / (2.0 * mu_brake * g)
    return s_reaction + s_brake

def required_TW_for_BFL_by_grid(ws, rho, CL_to, runway_length,
                                TW_min=0.01, TW_max=1.0, ngrid=400,
                                mu_roll=0.02, mu_brake=0.4, t_reaction=1.5,
                                airborne_factor=1.2, Vr_factor=1.2):
    """
    Grid-solve for T/W that makes (takeoff_to_obstacle - stop_distance) == runway_length.
    This is an approximate BFL formulation: we compute Vr = Vr_factor * V_stall,
    ground roll depends on acceleration (TW - mu_roll); airborne distance uses a simple scaling.
    Return numeric required T/W; np.nan if not found in bracket.
    """
    V_stall = stall_speed(ws, rho, CL_to)
    Vr = Vr_factor * V_stall
    s_stop = stop_distance_from_speed(Vr, mu_brake=mu_brake, t_reaction=t_reaction)

    TW_grid = np.linspace(TW_min, TW_max, ngrid)
    # compute takeoff distances for grid
    s_total = []
    for TW in TW_grid:
        s_ground = takeoff_ground_roll(Vr, TW, mu_roll=mu_roll)
        s_air = takeoff_airborne_distance(Vr, climb_angle_est=0.03, factor=airborne_factor)
        s_total.append(s_ground + s_air - s_stop)

    s_total = np.array(s_total)
    # target: s_total == runway_length
    diff = s_total - runway_length
    # find where sign changes or min abs
    if np.any(np.sign(diff[:-1]) != np.sign(diff[1:])):
        # locate first bracket and linear-interpolate for a better estimate
        idx = np.where(np.sign(diff[:-1]) != np.sign(diff[1:]))[0][0]
        x0, x1 = TW_grid[idx], TW_grid[idx+1]
        y0, y1 = diff[idx], diff[idx+1]
        # linear interpolation
        if (y1 - y0) != 0:
            TW_est = x0 - y0 * (x1 - x0) / (y1 - y0)
        else:
            TW_est = x0
        return float(TW_est)
    else:
        # no sign change -> return TW at minimal absolute difference (closest)
        idx = np.argmin(np.abs(diff))
        return float(TW_grid[idx])

# ----------------- Additional constraints -----------------
def required_TW_cruise(ws, rho, V, CD0, k):
    return drag_over_weight(ws, rho, V, CD0, k)

def required_TW_climb(ws, rho, V, CD0, k, climb_gradient):
    # simple: T/W = D/W + climb_gradient
    return required_TW_cruise(ws, rho, V, CD0, k) + climb_gradient

def required_TW_OEI(ws, rho, V, CD0, k, climb_gradient_req, n_engines=2, engines_failed=1):
    # D/W at condition, avail fraction after engine failure:
    D_W = drag_over_weight(ws, rho, V, CD0, k)
    avail_frac = (n_engines - engines_failed) / n_engines
    if avail_frac <= 0:
        return np.nan
    # T_total/W = (climb_grad + D/W) / avail_frac
    return (climb_gradient_req + D_W) / avail_frac

def required_TW_turn(ws, rho, V, CD0, k, load_factor=2.0):
    # sustained turn: CL scales with load factor, then D/W computed
    q = 0.5 * rho * V**2
    CL = 2.0 * ws * load_factor / (rho * V**2)
    CD = CD0 + k * CL**2
    return (q * CD) / ws

def required_TW_service_ceiling(ws, rho, V, CD0, k, ROC_req=0.5):
    # RoC = V * (T/W - D/W) -> T/W = D/W + ROC/V
    D_W = drag_over_weight(ws, rho, V, CD0, k)
    return D_W + ROC_req / V

# ----------------- Presets (SI units) -----------------
# To select a different aircraft: change ACTIVE_PRESET below (simple commenting/uncommenting).
presets = {
    "ATR72": {
        "MTOW_kg": 23000.0,
        "S_m2": 61.0,
        "CD0": 0.020,
        "k": 0.045,
        "CL_max_TO": 1.6,
        "CL_max_land": 2.0,
        "V_cruise_ms": 63.9,
        "alt_cruise_m": 6000,
        "runway_m": 1400.0,
        "climb_gradient": 0.024,
        "n_engines": 2
    },
    "Embraer E170": {
        "MTOW_kg": 37100.0,
        "S_m2": 54.46,
        "CD0": 0.017,
        "k": 0.045,
        "CL_max_TO": 1.7,
        "CL_max_land": 2.1,
        "V_cruise_ms": 230.0,
        "alt_cruise_m": 10000,
        "runway_m": 2000.0,
        "climb_gradient": 0.030,
        "n_engines": 2
    },
    "Airbus A320": {
        "MTOW_kg": 73500.0,
        "S_m2": 122.6,
        "CD0": 0.016,
        "k": 0.040,
        "CL_max_TO": 1.8,
        "CL_max_land": 2.2,
        "V_cruise_ms": 230.0,
        "alt_cruise_m": 11000,
        "runway_m": 2600.0,      # available runway length used for BFL calculation
        "climb_gradient": 0.030,
        "n_engines": 2
    },
    "Widebody (B777-like)": {
        "MTOW_kg": 300000.0,
        "S_m2": 427.8,
        "CD0": 0.015,
        "k": 0.035,
        "CL_max_TO": 2.0,
        "CL_max_land": 2.4,
        "V_cruise_ms": 250.0,
        "alt_cruise_m": 11000,
        "runway_m": 3000.0,
        "climb_gradient": 0.030,
        "n_engines": 2
    }
}

# ----------------- Choose active preset (A320 by default) -----------------
# To switch aircraft, set ACTIVE_PRESET to another key from 'presets' (comment/uncomment or edit).
ACTIVE_PRESET = "Airbus A320"
# Example for quick switching:
# ACTIVE_PRESET = "ATR72"
# ACTIVE_PRESET = "Embraer E170"
# ACTIVE_PRESET = "Widebody (B777-like)"

p = presets[ACTIVE_PRESET]

# Derived
W_N = p["MTOW_kg"] * g
S = p["S_m2"]
ws_actual = W_N / S            # N/m^2

# W/S sweep (SI); adjust limits as desired
ws_min = max(200.0, ws_actual*0.4)
ws_max = max(12000.0, ws_actual*2.5)
ws = np.linspace(ws_min, ws_max, 300)

# density values
rho_sea = rho_sl
rho_cruise = rho_sl * (0.7)    # approximate density at cruise for plotting convenience (you can compute ISA)

# Compute constraints for each ws
#  - cruise at altitude condition (approx using a chosen rho_cruise and V_cruise)
Vc = p["V_cruise_ms"]
CD0 = p["CD0"]
k = p["k"]
CL_to = p["CL_max_TO"]

tw_cruise = np.array([required_TW_cruise(w, rho_cruise, Vc, CD0, k) for w in ws])
tw_climb = np.array([required_TW_climb(w, rho_sea, Vc, CD0, k, p["climb_gradient"]) for w in ws])

# BFL: compute required T/W for available runway (grid solve)
runway_m = p["runway_m"]
tw_bfl = np.array([required_TW_for_BFL_by_grid(w, rho_sea, CL_to, runway_m,
                                               TW_min=0.01, TW_max=0.8, ngrid=600,
                                               mu_roll=0.02, mu_brake=0.4, t_reaction=1.5,
                                               airborne_factor=1.2, Vr_factor=1.2) for w in ws])

# OEI climb requirement (example target gradient 2.4% = 0.024)
tw_oei = np.array([required_TW_OEI(w, rho_sea, Vc, CD0, k, climb_gradient_req=0.024,
                                   n_engines=p.get("n_engines",2), engines_failed=1) for w in ws])

# Sustained turn (n=2)
tw_turn = np.array([required_TW_turn(w, rho_cruise, Vc, CD0, k, load_factor=2.0) for w in ws])

# Service ceiling requirement: RoC = 0.5 m/s
tw_service = np.array([required_TW_service_ceiling(w, rho_cruise, Vc, CD0, k, ROC_req=0.5) for w in ws])

# Landing W/S limit from approach speed and landing CL
V_app = p["V_cruise_ms"] * 0.32 + 65.0  # a rough approach speed heuristic or use p["V_app"] if defined
# Use explicit approach speed if you prefer; the preset can be extended to include "V_approach_ms".
CL_land = p["CL_max_land"] if "CL_max_land" in p else CL_to * 1.2
ws_landing_limit = 0.5 * rho_sea * CL_land * (V_app**2)

# Design point (actual)
tw_design = required_TW_cruise(ws_actual, rho_cruise, Vc, CD0, k)

# ----------------- Plot -----------------
fig, ax = plt.subplots()
ax.plot(ws, tw_cruise, label=f"{ACTIVE_PRESET} — cruise")
ax.plot(ws, tw_climb, linestyle='--', label=f"{ACTIVE_PRESET} — climb req")
ax.plot(ws, tw_bfl, linestyle='-.', label=f"{ACTIVE_PRESET} — BFL (runway {runway_m:.0f} m)")
ax.plot(ws, tw_oei, linestyle=':', label=f"{ACTIVE_PRESET} — OEI climb req (2.4%)")
ax.plot(ws, tw_turn, linestyle=(0, (1,5)), label=f"{ACTIVE_PRESET} — sustained turn n=2")
ax.plot(ws, tw_service, linestyle=(0, (3,1,1,1)), label=f"{ACTIVE_PRESET} — service ceiling (RoC 0.5 m/s)")

# Landing limit marker
# ax.axvline(ws_landing_limit, color='k', alpha=0.6, linewidth=1.0)
# ax.text(ws_landing_limit*1.02, ax.get_ylim()[1]*0.9, "Landing limit", rotation=90, va='top')

# Design point marker
ax.scatter([ws_actual], [tw_design], color='k', s=40, zorder=6)
ax.text(ws_actual*1.02, tw_design*1.02, f"Design point ({ACTIVE_PRESET})", fontsize=9)

# Axis labels (optionally convert W/S to imperial for display)
if use_imperial:
    # conversion 1 N/m^2 = 0.0208854342327 lb/ft^2
    conv = 0.0208854342327
    ax.set_xlabel("Wing loading W/S (lb/ft²)")
    # convert xtick labels numerically for readability
    xt = ax.get_xticks()
    ax.set_xticklabels([f"{x*conv:.1f}" for x in xt])
else:
    ax.set_xlabel("Wing loading W/S (N/m²)")

ax.set_ylabel("Thrust-to-weight ratio T/W (dimensionless)")
ax.set_title(f"CS-25 constraint diagram — {ACTIVE_PRESET} (approximate)")
ax.set_xlim(ws_min, ws_max)
ax.set_ylim(0.0, max(1.0, np.nanmax(np.hstack([tw_cruise, tw_climb, tw_bfl, tw_oei, tw_turn, tw_service]))*1.05))
ax.grid(True, which='both', linestyle=':')
ax.legend(fontsize=9, loc='upper right', ncol=1)

plt.tight_layout()
if save_png:
    outname = "constraint_a320.png"
    plt.savefig(outname, dpi=200)
    print(f"Saved diagram to {outname}")

plt.show()

# ----------------- Print a small numeric summary -----------------
print("\nPreset summary:")
print(f"  Aircraft: {ACTIVE_PRESET}")
print(f"  MTOW: {p['MTOW_kg']:.0f} kg, Wing area: {S:.1f} m^2, Design W/S: {ws_actual:.1f} N/m^2")
print(f"  BFL runway used: {runway_m:.0f} m")
print(f"  Landing W/S limit (approx): {ws_landing_limit:.1f} N/m^2")
if use_imperial:
    print("  (W/S axis displayed in lb/ft^2)")

# End of script
