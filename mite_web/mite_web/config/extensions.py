from flask_restx import Api

api = Api(
    title="MITE API",
    description="Rudimentary API to interact with MITE data repository.",
    doc="/api/",
    version="0.1",
)
