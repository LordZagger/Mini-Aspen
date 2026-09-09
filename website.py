from fastapi import FastAPI
from pydantic import BaseModel
import Zmodule260908
from fastapi.middleware.cors import CORSMiddleware
import io, sys

app = FastAPI(title="Mini-Aspen")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/docs/washing")
def docs_washing():
    return {"doc": Zmodule260908.washing_mccabe_thiele.__doc__}

@app.get("/docs/abs_strip_single")
def docs_abs_strip_single():
    return {"doc": Zmodule260908.absorption_stripping_single_stage.__doc__}

@app.get("/docs/abs_strip_multi")
def docs_abs_strip_multi():
    return {"doc": Zmodule260908.absorption_stripping_multi_stage.__doc__}

@app.get("/docs/flash")
def docs_flash():
    return {"doc": Zmodule260908.binary_flash_drum_sizing.__doc__}

@app.get("/docs/distillation")
def docs_distillation():
    return {"doc": Zmodule260908.binary_distillation_mccabe_thiele.__doc__}

@app.get("/docs/column")
def docs_column():
    return {"doc": Zmodule260908.column_diameter.__doc__}

@app.get("/docs/fugk")
def docs_fugk():
    return {"doc": Zmodule260908.FUGK.__doc__}

# Utility to capture printed output
def capture_output(func, *args, **kwargs):
    buffer = io.StringIO()
    sys.stdout = buffer
    func(*args, **kwargs)
    sys.stdout = sys.__stdout__
    return buffer.getvalue()


# -----------------------------
# 1. Washing McCabe–Thiele
# -----------------------------
class WashingInput(BaseModel):
    U: float
    O: float
    x_out: float
    x_in: float
    yNplus1: float = 0

@app.post("/washing")
def washing_endpoint(data: WashingInput):
    output = capture_output(
        Zmodule260908.washing_mccabe_thiele,
        data.U, data.O, data.x_out, data.x_in, data.yNplus1
    )
    return {"output": output}


# -----------------------------
# 2. Absorption/Stripping (Single Stage)
# -----------------------------
class AbsStripSingleInput(BaseModel):
    L_in: float
    V_in: float
    x_in: float
    y_in: float
    P: float
    H: float
    dilute: bool
    AbsStrip: bool
    linear: bool = False

@app.post("/abs_strip_single")
def abs_strip_single_endpoint(data: AbsStripSingleInput):
    result = Zmodule260908.absorption_stripping_single_stage(
        data.L_in, data.V_in, data.x_in, data.y_in,
        data.P, data.H, data.dilute, data.AbsStrip, data.linear
    )
    return {"result": result}


# -----------------------------
# 3. Absorption/Stripping (Multi Stage)
# -----------------------------
class AbsStripMultiInput(BaseModel):
    L_in: float
    V_in: float
    x_in: float
    y_in: float
    P: float
    H: float
    dilute: bool
    AbsStrip: bool
    x_out: float | None = None
    y_out: float | None = None

@app.post("/abs_strip_multi")
def abs_strip_multi_endpoint(data: AbsStripMultiInput):
    output = capture_output(
        Zmodule260908.absorption_stripping_multi_stage,
        data.L_in, data.V_in, data.x_in, data.y_in,
        data.P, data.H, data.dilute, data.AbsStrip,
        data.x_out, data.y_out
    )
    return {"output": output}


# -----------------------------
# 4. Flash Drum Sizing
# -----------------------------
class FlashInput(BaseModel):
    p_a: float
    p_b: float
    MW_a: float
    MW_b: float
    L: float
    V: float
    x_a: float
    y_a: float
    P: float
    T: float

@app.post("/flash_drum")
def flash_drum_endpoint(data: FlashInput):
    D = Zmodule260908.binary_flash_drum_sizing(
        data.p_a, data.p_b, data.MW_a, data.MW_b,
        data.L, data.V, data.x_a, data.y_a,
        data.P, data.T
    )
    return {"diameter_ft": D}


# -----------------------------
# 5. Binary Distillation (McCabe–Thiele)
# -----------------------------
class DistillationInput(BaseModel):
    xD: float
    xB: float
    a: float
    z: float
    full_mode: bool = True
    R: float | None = None
    F: float | None = None
    q: float | None = None
    D: float | None = None
    B: float | None = None
    L: float | None = None
    V: float | None = None
    Lp: float | None = None
    Vp: float | None = None
    factor: float | None = None

@app.post("/distillation")
def distillation_endpoint(data: DistillationInput):
    output = capture_output(
        Zmodule260908.binary_distillation_mccabe_thiele,
        data.xD, data.xB, data.a, data.z,
        data.full_mode, data.R, data.F, data.q,
        data.D, data.B, data.L, data.V,
        data.Lp, data.Vp, data.factor
    )
    return {"output": output}


# -----------------------------
# 6. Column Diameter
# -----------------------------
class ColumnDiameterInput(BaseModel):
    WL: float
    pL: float
    WV: float
    pV: float
    sigma: float
    Q: float
    f: float
    n: float
    spacing: int

@app.post("/column_diameter")
def column_diameter_endpoint(data: ColumnDiameterInput):
    d = Zmodule260908.column_diameter(
        data.WL, data.pL, data.WV, data.pV,
        data.sigma, data.Q, data.f, data.n, data.spacing
    )
    return {"diameter_ft": d}
