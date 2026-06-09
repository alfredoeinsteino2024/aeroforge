# ⚡ AeroForge

Browser-based aircraft aerodynamics simulator built with a C physics engine compiled to WebAssembly.

## Live Demo
Open `frontend/public/index.html` via a local server:
```cmd
cd frontend/public
npx serve .
```
Then open `http://localhost:3000`

## What it does
- Design any wing shape using parametric controls
- Real-time aerodynamic simulation (VLM solver)
- 3D pressure visualization with Three.js
- Animated streamlines showing airflow
- Stall warning at high angle of attack
- CL, CD, L/D ratio, lift force calculations

## Tech Stack
- **Physics engine:** C (Vortex Lattice Method solver)
- **WebAssembly:** Emscripten compiler
- **3D rendering:** Three.js
- **Frontend:** Vanilla HTML/CSS/JS

## Project Structure
aeroforge/
├── engine/
│   ├── src/vlm_solver.c     # C physics engine
│   ├── include/aero.h       # Data structures
│   └── build/               # Compiled WASM output
├── frontend/
│   └── public/index.html    # Main simulator UI
└── tests/
└── naca0012_baseline.c  # Physics validation

## Author
Toluwanimi Alfred — [@AlfredFadipe](https://x.com/AlfredFadipe)