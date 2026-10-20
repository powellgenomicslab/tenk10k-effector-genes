#!/usr/bin/env python3
"""Set a published Zenodo record's creators from zenodo/creators.json.

.zenodo.json, which Zenodo's GitHub integration reads for each code release, uses
Zenodo's legacy schema: one affiliation per creator, so it lists only the first.
After each GitHub release, run this on the new code record to restore every
author's affiliations (zenodo/creators.json, matching resources/authors.docx).
It opens an edit of the record, replaces the creators (same people, same order),
and saves it as a draft; with --publish it also publishes the edit. Editing
metadata does not create a new version.

    python3 zenodo/set_creators.py --record 23487290 [--publish]

Uses the same token file as upload_data.py (scope deposit:write; publishing
also needs deposit:actions). On NCI Gadi, run it on a login node or in copyq.
"""

import argparse
import json
import os
import sys
import urllib.error
import urllib.request

API = "https://zenodo.org/api"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from upload_data import read_token  # noqa: E402


def rdm(method, path, token, body=None):
    req = urllib.request.Request(f"{API}/{path}", method=method,
                                 data=None if body is None else json.dumps(body).encode())
    req.add_header("Accept", "application/vnd.inveniordm.v1+json")
    req.add_header("Content-Type", "application/json")
    req.add_header("Authorization", f"Bearer {token}")
    try:
        with urllib.request.urlopen(req) as resp:
            return json.loads(resp.read().decode() or "{}")
    except urllib.error.HTTPError as e:
        sys.exit(f"{method} {path} failed: HTTP {e.code}\n{e.read().decode()[:1500]}")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--record", required=True, help="record ID of the published record")
    p.add_argument("--token-file", default="~/.zenodo_token")
    p.add_argument("--publish", action="store_true", help="publish the edit")
    args = p.parse_args()

    token = read_token(args.token_file) or sys.exit(f"no token file {args.token_file}")
    creators = json.load(open(os.path.join(HERE, "creators.json")))

    try:  # open an edit draft; an already open one is reused
        rdm("POST", f"records/{args.record}/draft", token)
    except SystemExit:
        pass
    draft = rdm("GET", f"records/{args.record}/draft", token)
    meta = draft["metadata"]
    have = [c["person_or_org"]["name"] for c in meta["creators"]]
    want = [f'{c["person_or_org"]["family_name"]}, {c["person_or_org"]["given_name"]}' for c in creators]
    if have != want:
        sys.exit(f"creators differ from zenodo/creators.json; not changed\n record: {have}\n file:   {want}")
    meta["creators"] = creators
    # a licence is sent back by its id alone; Zenodo rejects the id with its display text
    meta["rights"] = [{"id": r["id"]} if r.get("id") else r for r in meta.get("rights", [])]
    body = {k: draft[k] for k in ("metadata", "access", "files", "pids", "custom_fields") if k in draft}
    saved = rdm("PUT", f"records/{args.record}/draft", token, body)
    if saved.get("errors"):
        sys.exit(f"not saved: {saved['errors']}")
    n = len({a.get("name") for c in saved["metadata"]["creators"] for a in c.get("affiliations", [])})
    print(f"record {args.record}: {len(creators)} creators, {n} distinct affiliations saved to the edit draft")
    if args.publish:
        rdm("POST", f"records/{args.record}/draft/actions/publish", token)
        print("published")
    else:
        print(f"review and publish at https://zenodo.org/uploads/{args.record}")


if __name__ == "__main__":
    main()
