"""Make the Dockerfile from `snakemake --containerize` reproducible.

Snakemake writes one `conda env create --file environment.yaml` per environment,
which re-solves the environment at build time and skips post-deploy scripts.
For every environment with a `<name>.linux-64.pin.txt` next to its YAML, this
rewrites the build to install the pinned explicit package list and run
`<name>.post-deploy.sh`.

Snakemake looks for each environment inside the container at /conda-envs/<hash>,
where <hash> covers the post-deploy script, the pin file for the platform it is
running on, and the YAML. The image is Linux, so the prefixes are recomputed
with the linux-64 pin; otherwise a Dockerfile generated on macOS names them by
the osx-arm64 pin and Snakemake on Linux cannot find them.

usage: python3 workflow/envs/pin_dockerfile.py Dockerfile
"""

import hashlib
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
    # same order as snakemake.deployment.conda.Env._get_hash
    md5 = hashlib.md5()
    if post.exists():
        md5.update(post.read_bytes())
    md5.update(pin.read_bytes())
    md5.update(yaml_path.read_bytes())
    text = text.replace(prefix, f"/conda-envs/{md5.hexdigest()}")

# link the GHCR package to its source repository (shown on the repository page
# and inheriting its access settings)
label = (
    "LABEL org.opencontainers.image.source="
    "https://github.com/powellgenomicslab/tenk10k-effector-genes\n"
)
if label not in text:
    text = text.replace('LABEL io.github.snakemake.containerized="true"\n',
                        'LABEL io.github.snakemake.containerized="true"\n' + label, 1)

path.write_text(text)
