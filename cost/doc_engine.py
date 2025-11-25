import json
import numpy as np
import pandas as pd
from datetime import datetime

# Global dictionary storing aircraft characteristics
AIRCRAFT = {
    "ATR72-600": {
        "mtow_kg": 23000, "seat_count": 70,
        "cruise_speed_kts": 275,
        "sfc_lb_hp_hr": 0.50, "installed_power_hp": 2750 * 2,
        "engine_count": 2, "list_price_musd": 26,
    },
    "Q400": {
        "mtow_kg": 29260, "seat_count": 78,
        "cruise_speed_kts": 360,
        "sfc_lb_hp_hr": 0.52, "installed_power_hp": 5071 * 2,
        "engine_count": 2, "list_price_musd": 34,
    },
    "E175": {
        "mtow_kg": 37200, "seat_count": 76,
        "cruise_speed_kts": 445,
        "sfc_lb_lbf_hr": 0.65, "installed_thrust_lbf": 14200 * 2,
        "engine_count": 2, "list_price_musd": 48,
    }
}

DEFAULT_FUEL_PRICE = 0.80


# ===============================================================
# LOAD CUSTOM AIRCRAFT FROM JSON
# ===============================================================
def load_custom_aircraft(json_content: str):
    data = json.loads(json_content)
    if not isinstance(data, dict):
        raise ValueError("JSON must contain dictionary of aircraft objects")

    for name, ac in data.items():
        validate_aircraft(name, ac)
        AIRCRAFT[name] = ac

    return list(data.keys())


# ===============================================================
# MANUAL AIRCRAFT ENTRY VALIDATION
# ===============================================================
def validate_aircraft(name, ac):
    required = ["mtow_kg", "seat_count", "cruise_speed_kts",
                "engine_count", "list_price_musd"]

    for r in required:
        if r not in ac:
            raise ValueError(f"{name}: Missing required field: {r}")

    turboprop = "sfc_lb_hp_hr" in ac and "installed_power_hp" in ac
    jet = "sfc_lb_lbf_hr" in ac and "installed_thrust_lbf" in ac

    if not (turboprop or jet):
        raise ValueError(
            f"{name} must define turboprop SFC+power OR jet SFC+thrust fields"
        )


# ===============================================================
# FUEL BURN MODEL
# ===============================================================
def fuel_burn_segmented(ac, mission_nm, fuel_price=DEFAULT_FUEL_PRICE):

    speed = ac["cruise_speed_kts"]
    cruise_hours = (mission_nm * 0.80) / speed
    climb_hours = 0.20
    descent_hours = 0.10

    if "sfc_lb_hp_hr" in ac:
        fuel_lb = ac["sfc_lb_hp_hr"] * ac["installed_power_hp"] * (climb_hours + cruise_hours + descent_hours)
    else:
        fuel_lb = ac["sfc_lb_lbf_hr"] * ac["installed_thrust_lbf"] * (climb_hours + cruise_hours + descent_hours)

    fuel_kg = fuel_lb * 0.453592
    return fuel_kg * fuel_price, fuel_kg


# ===============================================================
# COST MODELS
# ===============================================================
def maintenance_cost(ac, hours):
    return hours * (0.0009 * ac["mtow_kg"] + 75 * ac["engine_count"])

def crew_cost(hours):
    return 450 * hours

def nav_airport_cost(ac):
    mtow = ac["mtow_kg"]
    return 0.003 * mtow, 0.004 * mtow

def ownership_cost(ac, hours, utilization=2500, interest=0.08):
    price = ac["list_price_musd"] * 1e6
    return hours * ((price * interest) / utilization)


# ===============================================================
# MAIN DOC + CASM + RASM CALCULATOR
# ===============================================================
def compute_doc(ac_name, mission_nm, fuel_price=DEFAULT_FUEL_PRICE, avg_fare=120):

    ac = AIRCRAFT[ac_name]

    # Flight time
    block_hours = mission_nm / ac["cruise_speed_kts"] + 0.3

    fuel_cost, fuel_kg = fuel_burn_segmented(ac, mission_nm, fuel_price)
    maint = maintenance_cost(ac, block_hours)
    crew = crew_cost(block_hours)
    nav, apt = nav_airport_cost(ac)
    own = ownership_cost(ac, block_hours)

    components = {
        "Fuel": fuel_cost,
        "Maintenance": maint,
        "Crew": crew,
        "Navigation": nav,
        "Airport": apt,
        "Ownership": own
    }

    total = sum(components.values())

    # CASM (Cost per Available Seat Mile)
    asm = ac["seat_count"] * mission_nm
    casm = (total / asm) * 100  # in cents

    # RASM (Revenue per Available Seat Mile)
    revenue = avg_fare * ac["seat_count"]
    rasm = (revenue / asm) * 100  # cents

    return {
        "components": components,
        "total": total,
        "fuel_kg": fuel_kg,
        "block_hours": block_hours,
        "casm": casm,
        "rasm": rasm,
        "revenue": revenue,
        "per_hour": total / block_hours,
        "per_trip": total,
        "per_seat": total / ac["seat_count"],
    }
