"""Local configuration. Credentials never enter a run or an export."""

import os
from dataclasses import dataclass, field
from pathlib import Path

from dotenv import dotenv_values


@dataclass(frozen=True)
class Config:
    email: str = ""
    api_key: str = field(default="", repr=False)
    data_dir: Path = field(
        default_factory=lambda: (
            Path.home() / "Library" / "Application Support" / "KBioSuite" / "Kblast" / "runs"
        )
    )
    attempts: int = 3
    timeout: tuple[int, int] = (15, 90)
    blast_timeout_seconds: int = 6 * 60 * 60

    @classmethod
    def discover(cls, data_dir: Path | None = None) -> Config:
        values = dict(dotenv_values(Path.home() / ".config" / "kblast" / "ncbi.env"))
        values.update(dotenv_values(Path.cwd() / ".env"))
        values.update(os.environ)
        key = values.get("NCBI_API_KEY") or ""
        key_file = values.get("NCBI_API_KEY_FILE")
        if not key:
            path = Path(key_file).expanduser() if key_file else Path.home() / ".ncbi" / "api_key"
            if path.is_file():
                key = path.read_text().strip()
                if not key or any(c.isspace() for c in key):
                    raise ValueError("The NCBI key file must contain only the API key.")
        return cls(
            email=values.get("NCBI_EMAIL") or "",
            api_key=key.strip(),
            data_dir=data_dir
            or Path.home() / "Library" / "Application Support" / "KBioSuite" / "Kblast" / "runs",
        )
