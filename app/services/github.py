import logging

from async_lru import alru_cache
from fastapi import Request
from fastapi.concurrency import run_in_threadpool
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


@alru_cache(ttl=60)
async def get_kanban_cached(gh: Github):
    def fetch():
        repo = gh.get_repo(settings.repo_name)
        pulls = repo.get_pulls(state="open")
        return process_pulls(pulls)

    return await run_in_threadpool(fetch)


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


def create_pr(gh: Github, uuid: str):
    # TODO: complete implementation
    pass


def push_data(gh: Github, uuid: str, data: dict, name: str):
    # TODO: complete implementation
    pass
