import os
from pathlib import Path

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles

from app.api.utils import load_json

DATA_DIR = Path(os.environ.get("DATA_DIR", "/app/data"))

app = FastAPI()
app.mount("/data", StaticFiles(directory=DATA_DIR), name="data")


@app.get("/")
def read_root():
    version = load_json(DATA_DIR / "version_mite_data.json")

    return {"mite_data": version["version_mite_data"]}
