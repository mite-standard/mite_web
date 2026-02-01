from app.auth.passwords import verify_pw


def test_verify_pw_valid():
    assert verify_pw(
        plain="abc",
        hashed="$2b$12$inYu3eSuij5Xt9O8Pm75veEGPP/X.hlXpRNpu78YHSYGi8QHwzg.m",
    )


def test_verify_pw_invalid():
    assert not verify_pw(
        plain="abc",
        hashed="$2b$12$55OTrXkI7foWFPrD83wTh.7VQnCdW6Jga/D2sjKWhTqgfNVJqwfNC",
    )
