import bcrypt


def verify_pw(plain: str, hashed: str) -> bool:
    """Verify provided password with bcrypt-hashed one"""
    return bcrypt.checkpw(
        password=plain.encode("utf-8"), hashed_password=hashed.encode("utf-8")
    )
