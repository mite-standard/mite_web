import base64
import json
import sys
from pathlib import Path

import bcrypt
import pandas as pd


def main(csv: str):
    """Converts a csv of orcid-password pairs into a base-64-encoded string

    Input csv needs columns orcid,password

    """
    path = Path(csv)
    if not path.exists():
        raise FileNotFoundError

    df = pd.read_csv(path)
    required = {"orcid", "password"}
    if not required.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {required}")

    hashed = {}
    for _, row in df.iterrows():
        pw_hash = bcrypt.hashpw(
            row["password"].encode("utf-8"), bcrypt.gensalt(rounds=12)
        )

        hashed[row["orcid"]] = pw_hash.decode("utf-8")

    s = json.dumps(hashed)
    print(base64.b64encode(s.encode("utf-8")).decode("ascii"))


if __name__ == "__main__":
    main(csv=sys.argv[1])
