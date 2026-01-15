import logging

from fastapi import Request
from github import Auth, Github, PullRequest

from app.core.config import settings

logger = logging.getLogger(__name__)


def authenticate_pat() -> Github | None:
    """Authenticate with personal access"""
    if settings.env == "development":
        logger.warning("GitHub authentication in mode 'development' disabled")
        return

    auth = Auth.Token(settings.gh_token)
    return Github(auth=auth)


def get_github(request: Request) -> Github | None:
    if request.app.state.gh:
        return request.app.state.gh


def fake_pull() -> dict:
    # TODO: only for testing - remove later
    return {
        "Draft": [{"title": "Dummy", "url": "dummy", "created_at": "date"}],
        "In Review": [],
        "Reviewed": [],
    }


def process_pulls(pulls: list[PullRequest]) -> dict:
    """Process fetched pulls for submission page kanban board"""
    board = {"Draft": [], "In Review": [], "Reviewed": []}

    for pr in pulls:
        labels = {label.name for label in pr.labels}

        if pr.draft:
            column = "Draft"
        elif "reviewed" in labels:
            column = "Reviewed"
        else:
            column = "In Review"

        board[column].append(
            {"title": pr.title, "url": pr.html_url, "created_at": pr.created_at}
        )

    return board
