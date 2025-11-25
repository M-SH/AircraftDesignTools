import plotly.graph_objects as go

# Example DOC components for illustration
components = {
    "Fuel": 3000,
    "Maintenance": 1200,
    "Crew": 800,
    "Navigation": 200,
    "Airport": 300,
    "Ownership": 1500
}

# Create Pie Chart
fig = go.Figure(go.Pie(
    labels=list(components.keys()),
    values=list(components.values()),
    hole=0.4,
    hoverinfo="label+percent+value"
))

fig.update_layout(
    title_text="DOC Composition (Example Aircraft)",
    annotations=[dict(text="DOC", x=0.5, y=0.5, font_size=20, showarrow=False)]
)

# Add CASM/RASM as text annotations
casm = 0.145  # Example ¢/ASM
rasm = 0.180  # Example ¢/ASM

fig.add_annotation(
    x=0.5, y=-0.15,
    text=f"CASM: {casm:.2f} ¢/ASM<br>RASM: {rasm:.2f} ¢/ASM",
    showarrow=False,
    font=dict(size=14)
)

fig.show()
