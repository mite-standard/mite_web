import pytest

from app.core.config import Settings


def test_reviewers_valid():
    settings = Settings(
        reviewers="eyIwMDAwLTAwMDEtNjUzNC02NjA5IjogeyJwd19oYXNoZWQiOiAiJDJiJDEyJGhqRFlNZy9XUEN2V3RBYXFrNGxwZC5CUjMzRjk4cDlTZ05kZHYvZi5Bb1dlUDUzRmZPRUdPIiwgImdoX3RhZyI6ICJAZHVtbXkifX0="
    )
    assert settings.reviewer_gh_tags == ["@dummy"]
