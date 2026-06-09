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
- Design any wing shape using parametric controls (span, chord, sweep, dihedral, taper)
- Real-time aerodynamic simulation powered by a Vortex Lattice Method solver written in C
- 3D pressure visualization with Three.js — red = high pressure, blue = low pressure
- Animated streamlines showing airflow over the wing surface
- Stall warning at high angle of attack (AoA > 13°)
- Live CL, CD, L/D ratio, pitch moment, and actual lift force in Newtons

## Tech Stack
- **Physics engine:** C (Vortex Lattice Method solver compiled to WebAssembly)
- **WebAssembly compiler:** Emscripten 3.1.45
- **3D rendering:** Three.js r128
- **Frontend:** Vanilla HTML + CSS + JavaScript
- **Local server:** Node.js / npx serve

## Project Structure
aeroforge/
├── engine/
│   ├── src/
│   │   └── vlm_solver.c        # C physics engine — VLM + Biot-Savart solver
│   ├── include/
│   │   └── aero.h              # Data structures (WingGeometry, AeroResult, Panel)
│   └── build/
│       ├── aero.js             # Emscripten JS glue
│       └── aero.wasm           # Compiled WebAssembly binary
├── frontend/
│   ├── public/
│   │   ├── index.html          # Main simulator UI
│   │   ├── aero.js             # WASM glue (copied from engine/build)
│   │   └── aero.wasm           # WASM binary (copied from engine/build)
│   └── package.json
├── tests/
│   └── naca0012_baseline.c     # Physics validation against known NACA data
└── README.md

## How to recompile the physics engine
Requires [Emscripten](https://emscripten.org/docs/getting_started/downloads.html)

```cmd
cd C:\emsdk
emsdk_env.bat
cd C:\aeroforge
emcc engine\src\vlm_solver.c -I engine\include -O2 -s WASM=1 -s EXPORTED_FUNCTIONS="['_run_vlm_solver','_generate_mesh','_compute_lift','_compute_drag','_malloc','_free']" -s EXPORTED_RUNTIME_METHODS="['ccall','cwrap']" -s MODULARIZE=1 -s EXPORT_NAME="AeroEngine" -s ALLOW_MEMORY_GROWTH=1 -o engine\build\aero.js
copy /Y engine\build\aero.js frontend\public\aero.js
copy /Y engine\build\aero.wasm frontend\public\aero.wasm
```

## Physics
The solver implements the **Vortex Lattice Method (VLM)** — the same class of aerodynamic analysis used in early aircraft design tools. Each wing panel carries a horseshoe vortex whose influence on all other panels is computed via the Biot-Savart law. The resulting linear system is solved by Gaussian elimination, and lift is integrated using the Kutta-Joukowski theorem.

## Author
Toluwanimi Alfred — [@AlfredFadipe](https://x.com/AlfredFadipe) — [LinkedIn](https://linkedin.com/in/toluwanimialfred)