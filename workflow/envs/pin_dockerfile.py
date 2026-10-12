"""Make the Dockerfile from `snakemake --containerize` reproducible.

Snakemake writes one `conda env create --file environment.yaml` per environment,
which re-solves the environment at build time and skips post-deploy scripts.
For every environment with a `<name>.linux-64.pin.txt` next to its YAML, this
rewrites the build to install the pinned explicit package list and run
`<name>.post-deploy.sh`. The /conda-envs/<hash> prefixes are left unchanged,
because that is where Snakemake looks for them inside the container.

usage: python3 workflow/envs/pin_dockerfile.py Dockerfile
"""

import re
import sys
from pathlib import Path

path = Path(sys.argv[1])
text = path.read_text()

copies = re.findall(r"^COPY (\S+\.yaml) (/conda-envs/\w+)/environment\.yaml$", text, flags=re.M)
for yaml, prefix in copies:
    yaml_path = Path(yaml)
    pin = yaml_path.with_name(f"{yaml_path.stem}.linux-64.pin.txt")
    if not pin.exists():
        continue
    post = yaml_path.with_name(f"{yaml_path.stem}.post-deploy.sh")
    extra = f"COPY {pin} {prefix}/pin.txt\n"
    build = f"conda create --yes --prefix {prefix} --file {prefix}/pin.txt"
    if post.exists():
        extra += f"COPY {post} {prefix}/post-deploy.sh\n"
        build += f" && \\\n    conda run --prefix {prefix} bash {prefix}/post-deploy.sh"
    text = text.replace(
        f"COPY {yaml} {prefix}/environment.yaml\n",
        f"COPY {yaml} {prefix}/environment.yaml\n{extra}",
    )
    text = text.replace(
        f"conda env create --prefix {prefix} --file {prefix}/environment.yaml",
        build,
    )

path.write_text(text)
