from flask_restx import Api
from flask_sqlalchemy import SQLAlchemy

api = Api(
    title="MITE API",
    description="Rudimentary API to interact with MITE data repository.",
    doc="/api/v1",
    version="1",
)

db = SQLAlchemy()
