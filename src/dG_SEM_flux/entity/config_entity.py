from dataclass import dataclass
from pathlib import Path

@dataclass(frozen=True)
class dGElasticConfig:
    root_dir: Path
    node: str
    poly_deg: int
    time_integrator: str
