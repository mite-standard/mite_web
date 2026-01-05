from typing import Union

from fastapi import FastAPI

app = FastAPI()


@app.get("/")
def read_root():
    return {"Hello": "Wörld"}


@app.get("/item/{item_id}")
def item(item_id: int, q: str | None):
    """Returns item information"""
    return {"name": "Banana", "item_id": item_id, "opt": q}
