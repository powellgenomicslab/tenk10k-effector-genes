#!/usr/bin/env python3
"""Upload the data release to its reserved Zenodo draft and attach it to a community.

Uploads the files of the data record (data/zenodo/ by default) to an existing,
unpublished Zenodo draft, checking each file's MD5 against Zenodo's, and skips
files already uploaded unchanged. With --metadata it also sets the draft's
metadata from zenodo/data.json, keeping the reserved DOI. It then attaches the
draft to a Zenodo community (default: powellgenomicslab), so that it is
submitted there for review. It never publishes: check the draft in the browser
and click "Submit for review"; the record is published when a community
curator accepts it.

The token is read from a file (default ~/.zenodo_token, mode 600), never from
the command line, and is not printed. Create it at zenodo.org under
Account -> Applications -> Personal access tokens, with the deposit:write scope.

    python3 zenodo/upload_data.py --record 23406585 --dry-run
    python3 zenodo/upload_data.py --record 23406585 [--metadata]

On NCI Gadi, run it in a copyq job (the only queue with internet access).
"""

import argparse
import hashlib
import json
import os
import stat
import subprocess
import sys
import urllib.error
import urllib.request

API = "https://zenodo.org/api"
FILES = [
    "tenk10k_figure_inputs.tar.gz",
    "tenk10k_phase1_effector_genes.v5.parquet",
    "tenk10k_phase1_effector_genes.v5.md",
    "INVENTORY.tsv",
    "README.md",
]
METADATA_KEYS = ["title", "upload_type", "description", "creators", "license",
                 "keywords", "related_identifiers"]


def read_token(path):
    path = os.path.expanduser(path)
    if not os.path.exists(path):
        return None
    if os.stat(path).st_mode & (stat.S_IRWXG | stat.S_IRWXO):
        sys.exit(f"{path} is readable by others; run: chmod 600 {path}")
    return open(path).read().strip()


def api(method, url, token=None, body=None):
    url = url if url.startswith("http") else f"{API}/{url}"
    data = None if body is None else json.dumps(body).encode()
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("Accept", "application/json")
    if data is not None:
        req.add_header("Content-Type", "application/json")
    if token:
        req.add_header("Authorization", f"Bearer {token}")
    try:
        with urllib.request.urlopen(req) as resp:
            text = resp.read().decode()
            return json.loads(text) if text else {}
    except urllib.error.HTTPError as e:
        sys.exit(f"{method} {url.split('?')[0]} failed: HTTP {e.code}\n{e.read().decode()[:2000]}")


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def upload(path, url, token):
    # curl streams the file and shows progress; the token goes in a header read
    # from stdin, so it does not appear in the process list
    cmd = ["curl", "--silent", "--show-error", "--retry", "3", "-H", "@-",
           "--upload-file", path, "--write-out", "\n%{http_code}", url]
    res = subprocess.run(cmd, input=f"Authorization: Bearer {token}\n", universal_newlines=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    body, _, code = res.stdout.rpartition("\n")
    if res.returncode != 0 or not code.startswith("2"):
        sys.exit(f"upload of {os.path.basename(path)} failed (HTTP {code or '-'}):\n{body}{res.stderr}")
    return json.loads(body)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--record", required=True, help="record ID of the reserved draft, e.g. 23406585")
    p.add_argument("--files-dir", default="data/zenodo")
    p.add_argument("--token-file", default="~/.zenodo_token")
    p.add_argument("--community", default="powellgenomicslab", help="community slug ('' for none)")
    p.add_argument("--metadata", action="store_true", help="also set metadata from zenodo/data.json")
    p.add_argument("--dry-run", action="store_true", help="check everything, change nothing")
    args = p.parse_args()

    paths = [os.path.join(args.files_dir, f) for f in FILES]
    missing = [f for f in paths if not os.path.exists(f)]
    if missing:
        sys.exit(f"missing: {missing}")

    token = read_token(args.token_file)
    if token is None:
        sys.exit(f"no token file {args.token_file}")

    # 1. the draft
    dep = api("GET", f"deposit/depositions/{args.record}", token)
    if dep.get("submitted") or dep.get("state") == "done":
        sys.exit(f"record {args.record} is already published; upload a new version instead")
    doi = dep.get("metadata", {}).get("prereserve_doi", {}).get("doi")
    print(f"draft {args.record}: '{dep.get('title') or '(untitled)'}', reserved DOI {doi}")
    bucket = dep["links"]["bucket"]
    present = {o["key"]: o for o in api("GET", bucket, token).get("contents", [])}

    # 2. files
    for path in paths:
        name = os.path.basename(path)
        local = md5sum(path)
        remote = present.get(name, {}).get("checksum", "")
        size = os.path.getsize(path) / 1e6
        if remote == f"md5:{local}":
            print(f"  {name} ({size:,.1f} MB): already uploaded")
            continue
        if args.dry_run:
            print(f"  {name} ({size:,.1f} MB): would upload" + (" (replacing a different copy)" if remote else ""))
            continue
        print(f"  {name} ({size:,.1f} MB): uploading ...", flush=True)
        result = upload(path, f"{bucket}/{name}", token)
        if result.get("checksum") != f"md5:{local}":
            sys.exit(f"  {name}: checksum mismatch after upload ({result.get('checksum')} vs md5:{local})")
        print(f"  {name}: uploaded, MD5 verified")
    extra = sorted(set(present) - set(FILES))
    if extra:
        print(f"  note: the draft also holds {extra}; delete them in the browser if unintended")

    # 3. metadata (optional; keeps the reserved DOI and anything not in data.json)
    if args.metadata:
        meta = dict(dep.get("metadata", {}))
        new = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data.json")))
        meta.update({k: new[k] for k in METADATA_KEYS if k in new})
        if args.dry_run:
            print(f"metadata: would set {', '.join(k for k in METADATA_KEYS if k in new)} from data.json")
        else:
            api("PUT", f"deposit/depositions/{args.record}", token, {"metadata": meta})
            print("metadata: set from data.json")

    # 4. community
    if args.community:
        comm = api("GET", f"communities/{args.community}")
        title = comm["metadata"]["title"]
        if args.dry_run:
            print(f"community: would submit to '{title}' ({args.community})")
        else:
            api("PUT", f"records/{args.record}/draft/review", token,
                {"receiver": {"community": comm["id"]}, "type": "community-submission"})
            print(f"community: draft attached to '{title}' ({args.community})")

    print(f"\nNot published. Check the draft at https://zenodo.org/uploads/{args.record}, then click "
          f"'Submit for review'" + ("; the record is published when a curator accepts it." if args.community
                                     else "/'Publish'."))


if __name__ == "__main__":
    main()
