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


def is_genpept(acc: str, timeout: float = 7.0) -> bool:
    """Check if string is a valid genbank entry

    Returns:
        Bool if found/not found
        None if exception, e.g. timeout

    """
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    try:
        response = requests.get(
            url,
            params={
                "db": "protein",
                "term": f"{acc}[Accession]",
                "retmode": "json",
                "retmax": 1,
            },
            timeout=timeout,
        )
        response.raise_for_status()

        return int(response.json()["esearchresult"]["count"]) > 0

    except requests.Timeout:
        return False
    except requests.RequestException:
        return False
    except (ValueError, KeyError, TypeError):
        return False


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


def validate_nr_items(container: list, msg: str):
    """Raises ValueError if n=0 or invalid"""
    if len(container) == 0 or all(i == "" for i in container):
        raise ValueError(msg)
