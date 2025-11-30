import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

# -----------------------------
# Configuration
# -----------------------------
use_imperial = False  # Toggle for Imperial units (True = Imperial, False = SI)
allow_custom_input = True  # Enable interactive input
pio.renderers.default = "browser"

# Conversion factors
N_to_lbf = 0.224809
m2_to_ft2 = 10.7639
kg_to_lb = 2.20462

# Preset aircraft data (SI units)
presets = {
    "ATR72-600": {"weight": 22800 * 9.81, "wing_area": 61, "thrust": 2 * 2475 * 4.44822},
    "CRJ900": {"weight": 38330 * 9.81, "wing_area": 71, "thrust": 2 * 59600 * N_to_lbf},
    "A320": {"weight": 77000 * 9.81, "wing_area": 122.6, "thrust": 2 * 120000 * N_to_lbf},
    "B787-9": {"weight": 254700 * 9.81, "wing_area": 377, "thrust": 2 * 320000 * N_to_lbf}
}

# -----------------------------
# Aircraft selection
# -----------------------------
selected_aircraft = "ATR72-600"  # Default preset

if allow_custom_input:
    print("Available presets:")
    for name in presets.keys():
        print(f"- {name}")
    choice = input("Enter preset name or 'custom': ")
    if choice in presets:
        W0, S, T = presets[choice]["weight"], presets[choice]["wing_area"], presets[choice]["thrust"]
        selected_aircraft = choice
    else:
        try:
            if use_imperial:
                W0 = float(input("Enter aircraft weight (lbf): ")) / N_to_lbf
                S = float(input("Enter wing area (ft²): ")) / m2_to_ft2
                T = float(input("Enter thrust (lbf): ")) / N_to_lbf
            else:
                W0 = float(input("Enter aircraft weight (N): "))
                S = float(input("Enter wing area (m²): "))
                T = float(input("Enter thrust (N): "))
        except:
            W0, S, T = presets[selected_aircraft]["weight"], presets[selected_aircraft]["wing_area"], presets[selected_aircraft]["thrust"]
else:
    W0, S, T = presets[selected_aircraft]["weight"], presets[selected_aircraft]["wing_area"], presets[selected_aircraft]["thrust"]

# Derived parameters
wing_loading = W0 / S
thrust_to_weight = T / W0

# -----------------------------
# Generate design chart data
# -----------------------------
WL_range = np.linspace(1000, 8000, 200)

# CS-25 constraints (simplified models)
takeoff_TW = 0.00004 * WL_range + 0.2
landing_TW = 0.00003 * WL_range + 0.05
climb_TW = np.full_like(WL_range, 0.3)
cruise_TW = 0.00002 * WL_range + 0.1
service_ceiling_TW = np.full_like(WL_range, 0.25)
stall_speed_limit = 6000

# -----------------------------
# Plot with feasible region shading
# -----------------------------
fig = go.Figure()

# Add constraints with legend
fig.add_trace(go.Scatter(x=WL_range, y=takeoff_TW, mode='lines', name='Takeoff Constraint'))
fig.add_trace(go.Scatter(x=WL_range, y=landing_TW, mode='lines', name='Landing Constraint'))
fig.add_trace(go.Scatter(x=WL_range, y=climb_TW, mode='lines', name='Climb Constraint'))
fig.add_trace(go.Scatter(x=WL_range, y=cruise_TW, mode='lines', name='Cruise Constraint'))
fig.add_trace(go.Scatter(x=WL_range, y=service_ceiling_TW, mode='lines', name='Service Ceiling'))

# Stall speed limit line
fig.add_shape(type="line", x0=stall_speed_limit, y0=0, x1=stall_speed_limit, y1=0.6,
              line=dict(color="Purple", dash="dash"))

# Feasible region shading
upper_bound = np.minimum.reduce([takeoff_TW, landing_TW, climb_TW, cruise_TW, service_ceiling_TW])
fig.add_trace(go.Scatter(x=WL_range, y=upper_bound, fill='tozeroy', mode='none', name='Feasible Region', fillcolor='rgba(0,200,0,0.2)'))

# Design point
fig.add_trace(go.Scatter(x=[wing_loading], y=[thrust_to_weight], mode='markers',
                         name=f'{selected_aircraft} Design Point', marker=dict(size=10, color='red')))

# Axis labels
x_label = 'Wing Loading (lbf/ft²)' if use_imperial else 'Wing Loading (N/m²)'
fig.update_layout(title=f'CS-25 Conceptual Aircraft Design Chart - {selected_aircraft}',
                  xaxis_title=x_label,
                  yaxis_title='Thrust-to-Weight Ratio',
                  legend_title='Constraints',
                  width=1000, height=700)

fig.show()
# Save chart
#fig.write_image('cs25_aircraft_design_chart_presets.png')

# Print calculated performance metrics
print(f"Selected Aircraft: {selected_aircraft}")
print(f"Wing Loading: {wing_loading:.2f} N/m²")
print(f"Thrust-to-Weight Ratio: {thrust_to_weight:.3f}")
#("Chart saved as cs25_aircraft_design_chart_presets.png")
