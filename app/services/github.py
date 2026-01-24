import base64
import json
import logging
from typing import Union

from async_lru import alru_cache
from fastapi import HTTPException, Request
from fastapi.concurrency import run_in_threadpool
from github import Auth, Github, GithubException, PullRequest, Repository

from app.core.config import settings

logger = logging.getLogger(__name__)


def authenticate_pat() -> Github | None:
    """Authenticate with personal access"""
    if settings.env == "development":
        logger.warning("GitHub authentication in mode 'development' disabled")
        return

    auth = Auth.Token(settings.gh_token)
    return Github(auth=auth)


def get_github(request: Request) -> Union[Repository, None]:
    if request.app.state.repo:
        return request.app.state.repo


@alru_cache(ttl=60)
async def get_kanban_cached(repo: Repository):
    """Get overview of open PRs to visualize non-interactive Kanban"""

    def fetch():
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


async def create_pr(repo: Repository, branch: str, data: dict):
    """Create PR for new/modified entry"""

    def _push():
        try:
            try:
                repo.get_git_ref(f"heads/{branch}")
            except GithubException:
                base_ref = repo.get_git_ref("heads/main")
                repo.create_git_ref(
                    ref=f"refs/heads/{branch}",
                    sha=base_ref.object.sha,
                )

            repo.create_file(
                path=f"mite_data/data/{branch}.json",
                message=f"Add submission",
                content=json.dumps(data, indent=4),
                branch=branch,
            )

            repo.create_pull(
                title=f"Contributor submission {data['enzyme']['name']}",
                body="Automated draft submission - description TBA",
                head=branch,
                base="main",
                draft=True,
            )
        except Exception as e:
            raise HTTPException(
                400, detail=f"Error creating pull request: {e!s}"
            ) from e

    return await run_in_threadpool(_push)


async def draft_to_full(repo: Repository, branch: str) -> PullRequest:
    """Convert draft to full pull request"""
    body = f"""
A submission was performed via the MITE web portal and needs reviewing.

Submission ID: {branch}

## Review requested

{", ".join(settings.reviewer_gh_tags)}

## TODO Reviewers

- Review the entry [HERE](https://mite.bioinformatics.nl/submission/review/{branch}) (this link will only work for the current MITE release)

*This action was performed by `mite-bot`*
"""

    def _push():
        try:
            pr = None
            for p in repo.get_pulls(state="open"):
                if p.head.ref == branch:
                    pr = p
                    break

            if not pr:
                raise HTTPException(
                    status_code=404, detail=f"No open PR found for branch {branch}"
                )

            if not pr.draft:
                return pr

            pr.edit(body=body)

            if pr.draft:
                pr.mark_ready_for_review()

            return pr

        except GithubException as e:
            raise HTTPException(
                status_code=400,
                detail=f"GitHub error while updating PR: {e.data}",
            ) from e

    return await run_in_threadpool(_push)


async def upsert_json_file(repo: Repository, branch: str, data: dict, name: str):
    """Create or update file if existing"""

    def _push():
        path = f"mite_data/data/{name}.json"
        content = json.dumps(data, indent=4)

        try:
            existing = repo.get_contents(path, ref=branch)
            repo.update_file(
                path=path,
                message=f"Update submission",
                content=content,
                sha=existing.sha,
                branch=branch,
            )

        except GithubException as e:
            if e.status == 404:
                repo.create_file(
                    path=path,
                    message=f"Add submission",
                    content=content,
                    branch=branch,
                )
            else:
                raise

    return await run_in_threadpool(_push)


async def get_data(repo: Repository, branch: str) -> dict:
    """Pull data from remote"""

    def _fetch():
        try:
            contents = repo.get_contents(
                ref=branch, path=f"mite_data/data/{branch}.json"
            )
        except GithubException as e:
            raise HTTPException(404, detail=f"File does not exist: {e!s}") from e

        decoded = base64.b64decode(contents.content).decode("utf-8")
        return json.loads(decoded)

    return await run_in_threadpool(_fetch)


async def delete_file(repo: Repository, branch: str):
    """Pull data from remote"""

    def _push():
        path = f"mite_data/data/{branch}.json"
        try:
            file = repo.get_contents(path, ref=branch)
            repo.delete_file(
                path=path,
                message=f"Delete {path}",
                sha=file.sha,
                branch=branch,
            )

        except GithubException as e:
            if e.status == 404:
                return
            raise

    return await run_in_threadpool(_push)


async def add_pr_label(repo: Repository, branch: str, label: str):
    """Add 'reviewed' tag for merging on GitHub, optionally remove draft"""

    def _push():
        try:
            pr = None
            for p in repo.get_pulls(state="open"):
                if p.head.ref == branch:
                    pr = p
                    break

            if not pr:
                raise HTTPException(
                    status_code=404, detail=f"No open PR found for branch {branch}"
                )

            pr.add_to_labels(label)

            if pr.draft:
                pr.edit(draft=False)

        except GithubException as e:
            raise HTTPException(
                status_code=400,
                detail=f"GitHub error while adding PR label: {e.data}",
            ) from e

    return await run_in_threadpool(_push)
