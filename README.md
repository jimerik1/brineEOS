# Brine EOS – Density Calculation API

Flask + Gunicorn micro-service that implements the model in  
**SPE/IADC 16079, “Density Modeling for Pure and Mixed-Salt Brines as a Function of Composition, Temperature, and Pressure.”**

It calculates density for  

* Single-salt brines *(NaCl, KCl, CaCl₂, CaBr₂, ZnBr₂, ZnCl₂)*  
* Arbitrary mixed-salt brines  
* Pure water (CoolProp reference)  

across wide **P‑T** ranges and now supports both **MPa / K** *and* **bar / °C** input.

---

## ✨ Features
* Single & mixed brine solvers based on Brønsted–Guggenheim + Debye–Hückel
* CoolProp-based water density endpoint
* Strict JSON schema & validation
* Docker‑first deployment (`docker compose up --build`)
* `/healthz` endpoint for liveness probes
* Optional `pressure_unit` and `temperature_unit` fields (`MPa⇆bar`, `K⇆C`)
* Rotating file + stdout logging; optimised byte‑code (`PYTHONOPTIMIZE=1`)

---

## 🛠 Quick start

```bash
# clone
git clone https://github.com/<you>/brineEOS.git
cd brineEOS

# build & run (will expose 5099 on the host)
docker compose up --build
```

Logs:

```
Brine Density API master ready ‣ http://0.0.0.0:5099 → container :5000
```

Health probe:

```bash
curl http://localhost:5099/healthz
# → {"status":"ok"}
```

---

## 🔌 Endpoints

| Verb | Path | Purpose |
|------|------|---------|
| **POST** | `/api/v1/calculate_density` | Single‑salt or mixed‑salt brine |
| **POST** | `/api/v1/calculate_water_density` | Pure‑water density (CoolProp) |
| **GET**  | `/healthz` | Liveness check |

### 1  Request body (brine)

```jsonc
{
  "brine_type" : "single" | "mixed" | "NaCl" | "KCl" | "CaCl2" | "CaBr2" | "ZnBr2" | "ZnCl2",

  "salt_composition" : {           // required for mixed, optional for single
    "CaCl2" : 24.0,
    "CaBr2" : 25.3
  },

  "base_density" : 1200.0,         // required for single if no composition

  "pressure_interval" : [0.1, 150.0],
  "pressure_resolution": 10.0,
  "pressure_unit" : "MPa",         // "MPa" (default) or "bar"

  "temperature_interval" : [298.15, 447.0],
  "temperature_resolution": 10.0,
  "temperature_unit" : "K"         // "K" (default) or "C"
}
```

### 2  Request body (water)

Same numeric fields, plus optional unit keys.

### 3  Units recap

| Field | Default | Also accepted |
|-------|---------|---------------|
| `pressure_unit` | **MPa** | `bar` |
| `temperature_unit` | **K** | `C` |

Values are converted internally; the solver always works in MPa & K.

### 4  Typical response

```jsonc
{
  "metadata": {
    "brine_type": "mixed",
    "pressure_points": [1.0, 11.0, …],
    "temperature_points": [295.0, 320.0, …],
    "units": { "density":"kg/m³", "pressure":"MPa", "temperature":"K" }
  },
  "densities": {
    "1.00": { "295.00": 1574.0, "320.00": 1559.2, … },
    "11.00": { … }
  }
}
```

---

## 🖥 Local development

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt

export FLASK_DEBUG=1
gunicorn -c gunicorn_conf.py --bind 0.0.0.0:5000 --workers 2 app:create_app()
```

Set `API_DEBUG=1` to enable in‑code debug prints.

---

## 🧪 Testing

```bash
pytest
```

CI runs `pytest`, `ruff`, `mypy`, and `pip-audit` on every PR.

---

## 🗄 Folder structure

```
brineEOS/
│
├── app.py                 # Flask factory
├── calculators/           # single, mixed, water calculators
├── config/                # salt parameter tables
├── routes/                # Flask blueprints
├── utils/                 # validators, converters, CoolProp helper
├── tests/                 # pytest suite
├── Dockerfile
└── docker-compose.yml
```

---

## 📜 License

MIT – see [LICENSE](LICENSE).

**Reference**  
N.P. Kemp & D.C. Thomas, “Density Modeling for Pure and Mixed‑Salt Brines as a Function of Composition, Temperature, and Pressure,” *SPE/IADC 16079*.

Questions or ideas? Open an issue or start a discussion.
