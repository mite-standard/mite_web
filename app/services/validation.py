import re
from pathlib import Path

import requests


def validate_submitter(v: str) -> None:
    """Validate submitter ID

    Args:
        v: ORCID or anonymous

    Raises:
        ValueError
    """
    pattern = re.compile(r"^\d{4}-\d{4}-\d{4}-\d{3}[\dX]$|^A{24}$")
    if not re.fullmatch(pattern, v):
        raise ValueError("Not a valid ORCID/not left blank.")


def check_exists(p: Path) -> None:
    """Check if file exists

    Args:
        p: the file path

    Raises:
        FileNotFoundError
    """
    if not p.exists():
        raise FileNotFoundError(f"MITE accession does not exist {p}.")


def check_mite_id(s: str) -> None:
    """Check if file exists

    Args:
        s: the mite accession ID

    Raises:
        ValueError
    """
    if not re.fullmatch(re.compile(r"^MITE(\d{7})$"), s):
        raise ValueError("MITE Accession invalid.")


def is_genpept(acc: str, timeout: float = 7.0) -> bool | None:
    """Check if string is a valid genbank entry

    Returns:
        Bool if found/not found
        None if exception, e.g. timeout

    """
    try:
        r = requests.head(
            url=f"https://www.ncbi.nlm.nih.gov/protein/{acc}",
            timeout=timeout,
            allow_redirects=False,
        )
        return r.status_code == 200
    except requests.RequestException:
        return None


def is_uniprot(acc: str, timeout: float = 3.0) -> bool | None:
    """Check if string is an existing uniprot entry"""
    try:
        if acc.startswith("UPI"):
            url = f"https://rest.uniprot.org/uniparc/{acc}.fasta"
        else:
            url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        r = requests.head(url=url, timeout=timeout, allow_redirects=True)
        return r.status_code == 200
    except requests.RequestException:
        return None
