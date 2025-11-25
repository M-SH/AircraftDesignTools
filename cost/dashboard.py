from dash import Dash, html, dcc, Input, Output, State, callback_context
import plotly.express as px
import numpy as np
import pandas as pd

from doc_engine import AIRCRAFT, compute_doc, load_custom_aircraft, validate_aircraft

# Fixed colors for DOC components
DOC_COLORS = {
    "Fuel": "#1f77b4",         # blue
    "Crew": "#ff7f0e",         # orange
    "Maintenance": "#2ca02c",  # green
    "Navigation": "#d62728",   # red
    "Airport": "#9467bd",      # purple
    "Ownership": "#8c564b"     # brown
}

app = Dash(__name__)
app.title = "Aircraft DOC Dashboard"

# --------------------------------------------------------------
# UI Layout
# --------------------------------------------------------------
app.layout = html.Div([
    html.H1("Aircraft DOC Dashboard"),

    # ----------------------------------------------------------
    # JSON Upload
    # ----------------------------------------------------------
    html.H3("Upload Custom Aircraft (JSON)"),
    dcc.Upload(
        id="upload_json",
        children=html.Div(["Drag and drop or click to upload JSON"]),
        style={
            "width": "40%",
            "height": "60px",
            "lineHeight": "60px",
            "borderWidth": "1px",
            "borderStyle": "dashed",
            "borderRadius": "5px",
            "textAlign": "center",
            "marginBottom": "10px"
        },
    ),
    html.Div(id="upload_status", style={"color": "green", "marginBottom": "20px"}),

    html.Hr(),

    # ----------------------------------------------------------
    # Manual Aircraft Entry Panel
    # ----------------------------------------------------------
    html.H3("Add / Modify Aircraft Manually"),

    html.Div([
        html.Div([
            html.Label("Aircraft Name"),
            dcc.Input(id="ac_name", type="text", placeholder="MyAircraft", style={"width": "90%"}),

            html.Label("MTOW (kg)"),
            dcc.Input(id="ac_mtow", type="number", style={"width": "90%"}),

            html.Label("Seats"),
            dcc.Input(id="ac_seats", type="number", style={"width": "90%"}),

            html.Label("Cruise Speed (kts)"),
            dcc.Input(id="ac_speed", type="number", style={"width": "90%"})
        ], style={"width": "25%", "display": "inline-block", "verticalAlign": "top"}),

        html.Div([
            html.Label("Engine Count"),
            dcc.Input(id="ac_engines", type="number", style={"width": "90%"}),

            html.Label("List Price (MUSD)"),
            dcc.Input(id="ac_price", type="number", style={"width": "90%"}),

            html.H4("Turboprop Parameters (optional)"),
            html.Label("SFC (lb/hp/hr)"),
            dcc.Input(id="ac_sfc_tp", type="number", style={"width": "90%"}),

            html.Label("Installed Power (hp)"),
            dcc.Input(id="ac_power", type="number", style={"width": "90%"})
        ], style={"width": "25%", "display": "inline-block", "verticalAlign": "top"}),

        html.Div([
            html.H4("Jet Parameters (optional)"),
            html.Label("SFC (lb/lbf/hr)"),
            dcc.Input(id="ac_sfc_jet", type="number", style={"width": "90%"}),

            html.Label("Installed Thrust (lbf)"),
            dcc.Input(id="ac_thrust", type="number", style={"width": "90%"}),

            html.Button("Add / Update Aircraft", id="add_ac_btn",
                        style={"marginTop": "20px", "width": "90%"}),

            html.Div(id="manual_status", style={"color": "blue", "marginTop": "10px"})
        ], style={"width": "25%", "display": "inline-block", "verticalAlign": "top"})
    ]),

    html.Hr(),

    # ----------------------------------------------------------
    # Main Inputs
    # ----------------------------------------------------------
    html.Label("Select Aircraft"),
    dcc.Dropdown(list(AIRCRAFT.keys()), list(AIRCRAFT.keys())[0], id="aircraft"),

    html.Br(),
    html.Label("Mission Length (nm)"),
    dcc.Slider(100, 1200, 50, value=300, id="mission_nm"),

    html.Br(),
    html.Label("Fuel Price (USD/kg)"),
    dcc.Slider(0.5, 2.0, 0.05, value=0.80, id="fuel_price"),

    html.Br(),
    html.Label("Average Fare (USD)"),
    dcc.Slider(50, 400, 10, value=120, id="avg_fare"),

    html.H2("Key Metrics"),
    html.Div(id="metrics_box", style={"fontSize": "18px", "marginBottom": "20px"}),

    dcc.Graph(id="pie_chart"),
    dcc.Graph(id="bar_chart"),
    dcc.Graph(id="heatmap")
])


# --------------------------------------------------------------
# Combined callback for both JSON uploads and manual aircraft entry
# --------------------------------------------------------------
@app.callback(
    Output("upload_status", "children"),
    Output("manual_status", "children"),
    Output("aircraft", "options"),
    Input("upload_json", "contents"),
    Input("add_ac_btn", "n_clicks"),
    State("ac_name", "value"),
    State("ac_mtow", "value"),
    State("ac_seats", "value"),
    State("ac_speed", "value"),
    State("ac_engines", "value"),
    State("ac_price", "value"),
    State("ac_sfc_tp", "value"),
    State("ac_power", "value"),
    State("ac_sfc_jet", "value"),
    State("ac_thrust", "value"),
    prevent_initial_call=True
)
def handle_aircraft_updates(
        uploaded_json,
        add_clicked,
        name, mtow, seats, speed, engines, price,
        sfc_tp, power, sfc_jet, thrust):

    triggered = callback_context.triggered[0]["prop_id"]

    upload_msg = ""
    manual_msg = ""

    # ----------------------------------------------------------
    # JSON UPLOAD HANDLING
    # ----------------------------------------------------------
    if "upload_json" in triggered and uploaded_json is not None:
        import base64
        _, content_string = uploaded_json.split(",")
        decoded = base64.b64decode(content_string).decode("utf-8")

        try:
            added = load_custom_aircraft(decoded)
            upload_msg = f"Loaded aircraft: {', '.join(added)}"
        except Exception as e:
            upload_msg = f"Error loading JSON: {e}"

    # ----------------------------------------------------------
    # MANUAL AIRCRAFT ENTRY HANDLING
    # ----------------------------------------------------------
    if "add_ac_btn" in triggered:
        if not name:
            manual_msg = "Error: Aircraft name is required"
        else:
            ac = {
                "mtow_kg": mtow,
                "seat_count": seats,
                "cruise_speed_kts": speed,
                "engine_count": engines,
                "list_price_musd": price
            }

            if sfc_tp and power:
                ac["sfc_lb_hp_hr"] = sfc_tp
                ac["installed_power_hp"] = power

            if sfc_jet and thrust:
                ac["sfc_lb_lbf_hr"] = sfc_jet
                ac["installed_thrust_lbf"] = thrust

            try:
                validate_aircraft(name, ac)
                AIRCRAFT[name] = ac
                manual_msg = f"Aircraft '{name}' added/updated successfully."
            except Exception as e:
                manual_msg = f"Error: {e}"

    return upload_msg, manual_msg, list(AIRCRAFT.keys())


# --------------------------------------------------------------
# Metrics Box (CASM, RASM, DOC)
# --------------------------------------------------------------
@app.callback(
    Output("metrics_box", "children"),
    Input("aircraft", "value"),
    Input("mission_nm", "value"),
    Input("fuel_price", "value"),
    Input("avg_fare", "value")
)
def update_metrics(ac, nm, fp, fare):
    r = compute_doc(ac, nm, fp, fare)

    return html.Div([
        html.Div(f"Total DOC: ${r['total']:.0f}"),
        html.Div(f"Block Time: {r['block_hours']:.2f} h"),
        html.Div(f"Fuel Burn: {r['fuel_kg']:.0f} kg"),
        html.Div(f"CASM: {r['casm']:.2f} ¢ per ASM"),
        html.Div(f"RASM: {r['rasm']:.2f} ¢ per ASM"),
        html.Div(f"Revenue per Flight: ${r['revenue']:.0f}")
    ])


# --------------------------------------------------------------
# Pie Chart — DOC composition
# --------------------------------------------------------------
@app.callback(
    Output("pie_chart", "figure"),
    Input("aircraft", "value"),
    Input("mission_nm", "value"),
    Input("fuel_price", "value"),
    Input("avg_fare", "value"))
def update_pie(ac, nm, fp, fare):
    r = compute_doc(ac, nm, fp, fare)
    fig = px.pie(
        names=list(r["components"].keys()),
        values=list(r["components"].values()),
        title=f"DOC Composition — {ac}",
        color=list(r["components"].keys()),
        color_discrete_map=DOC_COLORS
    )
    return fig



# --------------------------------------------------------------
# Bar Chart — compares all aircraft
# --------------------------------------------------------------
@app.callback(
    Output("bar_chart", "figure"),
    Input("mission_nm", "value"),
    Input("fuel_price", "value"),
    Input("avg_fare", "value"))
def update_bar(nm, fp, fare):
    rows = []

    for ac in AIRCRAFT.keys():
        result = compute_doc(ac, nm, fp, fare)
        for comp, value in result["components"].items():
            rows.append({
                "Aircraft": ac,
                "Component": comp,
                "Cost (USD)": value
            })

    df = pd.DataFrame(rows)
    fig = px.bar(
        df,
        x="Aircraft",
        y="Cost (USD)",
        color="Component",
        title="DOC Composition Comparison Across Aircraft",
        text_auto=True,
        color_discrete_map=DOC_COLORS
    )
    fig.update_layout(barmode='stack')
    return fig


# --------------------------------------------------------------
# Heatmap — Sensitivity to Fuel Price & Utilization
# --------------------------------------------------------------
@app.callback(
    Output("heatmap", "figure"),
    Input("aircraft", "value"),
    Input("mission_nm", "value"))
def update_heatmap(ac, nm):

    fp_range = np.linspace(0.5, 1.5, 10)
    util_range = np.linspace(1500, 3500, 10)

    z = []
    for fp in fp_range:
        row = []
        for util in util_range:
            doc = compute_doc(ac, nm, fp)
            own_adj = (AIRCRAFT[ac]["list_price_musd"] * 1e6 * 0.08) / util
            total_adj = doc["total"] - doc["components"]["Ownership"] + own_adj * doc["block_hours"]
            row.append(total_adj)
        z.append(row)

    return px.imshow(
        z,
        x=util_range,
        y=fp_range,
        labels={"x": "Utilization (hours/year)", "y": "Fuel price (USD/kg)"},
        title=f"DOC Sensitivity Heatmap — {ac}",
        color_continuous_scale="Viridis"
    )


# --------------------------------------------------------------
# Run App
# --------------------------------------------------------------
if __name__ == "__main__":
    app.run(debug=True)
