"""
STRING Values/Ranks Enrichment via REST API
https://string-db.org/help/api/#valuesranks-enrichment-api

Snakemake script for rank/value-based enrichment analysis using STRING's
asynchronous values/ranks enrichment API endpoints.

Workflow:
  1. Obtain / reuse a persistent API key  (cached in ~/.string_api_key)
  2. Read phenotype gene-set files  (two-column: gene<TAB>value, no header)
  3. Optionally pre-map identifiers to STRING IDs for faster processing
  4. Submit one enrichment job per phenotype file
  5. Poll all jobs until completion
  6. Download, combine and write results as a single TSV

The output schema is aligned with the gprofiler.R / stringdb.R scripts so
that aggregate_enrichment can consume it without changes.

Compatible with the enrichment_by_method_py rule in enrichment.smk.
"""

import io
import json
import logging
import os
import sys
import time
from pathlib import Path

import pandas as pd
import requests

# ── Snakemake interface ──────────────────────────────────────────────
INPUT     = snakemake.input
OUTPUT    = snakemake.output
PARAMS    = snakemake.params
LOG       = snakemake.log
WILDCARDS = snakemake.wildcards

# ── Logging ──────────────────────────────────────────────────────────
log_file = str(LOG[0])
os.makedirs(os.path.dirname(log_file), exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler(log_file, mode="w"),
        logging.StreamHandler(sys.stderr),
    ],
)
log = logging.getLogger(__name__)

# ── Parameters (with defaults matching enrichment.smk) ───────────────
species         = int(getattr(PARAMS, "species", 9606))
string_version  = str(getattr(PARAMS, "string_version", "12.0"))
min_genes       = int(getattr(PARAMS, "min_genes", 5))
fdr_threshold   = float(getattr(PARAMS, "fdr_threshold", 0.05))
caller_id       = str(getattr(PARAMS, "caller_identity", "tenk10k_smr"))
poll_interval   = int(getattr(PARAMS, "poll_interval", 30))
max_wait        = int(getattr(PARAMS, "max_wait", 7200))
api_key_path    = Path(os.path.expanduser(
    str(getattr(PARAMS, "api_key_file", ".string_api_key"))))
premap_ids      = bool(getattr(PARAMS, "premap_ids", False))
rank_direction  = int(getattr(PARAMS, "rank_direction", 0))
# rank_direction: 0 = both extremes (default), +1 = top only, -1 = bottom only

# Submission batching: pause after every N successful submissions to avoid
# flooding STRING's queue with many large jobs at once (which causes
# server-side 'globalenrichment.py error code 1' failures).
submit_batch_size  = int(getattr(PARAMS, "submit_batch_size", 20))
submit_batch_pause = int(getattr(PARAMS, "submit_batch_pause", 60))

STRING_API = "https://string-db.org/api"

log.info("Start: STRING Values/Ranks Enrichment API")
log.info(f"  dir_gene_set:   {INPUT.dir_gene_set}")
log.info(f"  gene_universe:  {INPUT.gene_universe}")
log.info(f"  output:         {OUTPUT.enrich}")
log.info(f"  species:        {species}")
log.info(f"  string_version: {string_version}")
log.info(f"  api_url:        {STRING_API}")
log.info(f"  min_genes:      {min_genes}")
log.info(f"  fdr_threshold:  {fdr_threshold}")
log.info(f"  poll_interval:  {poll_interval}s")
log.info(f"  max_wait:       {max_wait}s")
log.info(f"  premap_ids:     {premap_ids}")


# ══════════════════════════════════════════════════════════════════════
# API helpers
# ══════════════════════════════════════════════════════════════════════

def _request(method, endpoint, *, params=None, data=None, fmt="json",
             retries=3, timeout=300, politeness_delay=1.0,
             retry_on_client_error=True):
    """Fire a GET or POST to STRING API with retry and politeness delay.

    Parameters
    ----------
    politeness_delay : float
        Seconds to sleep before each attempt.  STRING asks for ≥1 s between
        heavy calls (submit, network, mapping).  Lightweight status-check
        endpoints can use a shorter delay (e.g. 0.2 s) to avoid spending
        minutes just polling when many jobs are in flight.
    retry_on_client_error : bool
        If False, 4xx errors are raised immediately without retrying.
        Useful for status-check endpoints where a 400 is not transient
        and retrying just wastes time.
    """
    url = f"{STRING_API}/{fmt}/{endpoint}"
    for attempt in range(retries):
        time.sleep(politeness_delay)
        try:
            if method == "GET":
                resp = requests.get(url, params=params, timeout=timeout)
            else:
                resp = requests.post(url, data=data, timeout=timeout)
            # If client error and we shouldn't retry, raise immediately
            if not retry_on_client_error and 400 <= resp.status_code < 500:
                resp.raise_for_status()
            resp.raise_for_status()
            return resp.json() if fmt == "json" else resp.text
        except requests.exceptions.HTTPError as exc:
            # Don't retry 4xx unless told to
            if not retry_on_client_error and exc.response is not None \
                    and 400 <= exc.response.status_code < 500:
                raise
            if attempt < retries - 1:
                # Use longer backoff for 5xx (server errors like 502)
                is_5xx = (exc.response is not None
                          and exc.response.status_code >= 500)
                wait = (10 * (attempt + 1)) if is_5xx else 2 ** (attempt + 1)
                log.warning(
                    f"  Retry {attempt + 1}/{retries} {endpoint}: {exc} "
                    f"(wait {wait}s)"
                )
                time.sleep(wait)
            else:
                raise
        except Exception as exc:
            if attempt < retries - 1:
                wait = 2 ** (attempt + 1)
                log.warning(
                    f"  Retry {attempt + 1}/{retries} {endpoint}: {exc} "
                    f"(wait {wait}s)"
                )
                time.sleep(wait)
            else:
                raise


def api_get(endpoint, params=None, **kw):
    return _request("GET", endpoint, params=params, **kw)


def api_post(endpoint, body, **kw):
    return _request("POST", endpoint, data=body, **kw)


# ── API key management ───────────────────────────────────────────────

def get_api_key():
    """Retrieve a cached STRING API key, or request and cache a new one.

    New keys activate within ~30 minutes.  The script retries automatically
    on submission if the key is not yet active.
    """
    if api_key_path.exists():
        key = api_key_path.read_text().strip()
        if key:
            log.info(f"  Using cached API key from {api_key_path}")
            return key

    log.info("  No API key found — requesting new key from STRING …")
    log.info(f"  TIP: to avoid activation delay, pre-generate a key:")
    log.info(f"    curl -s '{STRING_API}/json/get_api_key' | "
             f"python3 -c \"import sys,json; "
             f"print(json.load(sys.stdin)[0]['api_key'])\" "
             f"> {api_key_path}")

    result = api_get("get_api_key")
    key = result[0]["api_key"]
    api_key_path.parent.mkdir(parents=True, exist_ok=True)
    api_key_path.write_text(key + "\n")
    os.chmod(str(api_key_path), 0o600)
    log.info(f"  New API key saved to {api_key_path}")
    log.info("  NOTE: new keys take up to 30 min to activate; "
             "the script retries automatically when submitting jobs.")
    return key


# ── Identifier mapping (optional speed-up) ───────────────────────────

def map_identifiers(genes, batch_size=2000):
    """Map gene symbols → STRING IDs.  Returns {gene_name: stringId}."""
    mapping = {}
    for i in range(0, len(genes), batch_size):
        batch = genes[i : i + batch_size]
        try:
            txt = api_post(
                "get_string_ids",
                {
                    "identifiers": "\r".join(batch),
                    "species": species,
                    "limit": 1,
                    "echo_query": 1,
                    "caller_identity": caller_id,
                },
                fmt="tsv",
            )
            if txt.strip():
                df = pd.read_csv(io.StringIO(txt), sep="\t")
                for _, row in df.iterrows():
                    mapping[str(row["queryItem"])] = str(row["stringId"])
        except Exception as exc:
            log.warning(f"  map_identifiers batch {i}: {exc}")
    return mapping


# ── Read gene-set file ───────────────────────────────────────────────

def read_gene_file(path):
    """Read a gene-set file.

    Expects two TAB-separated columns: ``gene<TAB>value`` with NO header.
    Returns ``(identifiers_string, n_genes, is_ranked)``.

    * *identifiers_string*: newline-joined ``gene\\tvalue`` lines, ready
      for the ``identifiers`` parameter of the submit endpoint.
    * *is_ranked*: ``True`` when the file is valid two-column input.
    """
    text = Path(path).read_text()
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        return None, 0, False

    # Detect two-column (gene <TAB or space> numeric_value)
    sep = "\t" if "\t" in lines[0] else None
    parts = lines[0].split(sep, maxsplit=1)
    is_ranked = False
    if len(parts) >= 2:
        try:
            float(parts[1].split()[0])
            is_ranked = True
        except (ValueError, IndexError):
            pass

    if not is_ranked:
        return None, len(lines), False

    # Normalise every line to tab-separated
    normalised = []
    for ln in lines:
        fields = ln.split(sep, maxsplit=1)
        if len(fields) >= 2:
            normalised.append(f"{fields[0]}\t{fields[1]}")
    return "\n".join(normalised), len(normalised), True


# ── Submit a job ─────────────────────────────────────────────────────

def submit_job(identifiers_str, api_key, max_auth_wait=1800):
    """Submit a values/ranks enrichment job.

    Handles newly created API keys that are not yet activated by
    retrying for up to *max_auth_wait* seconds.
    """
    t0 = time.time()
    while True:
        result = api_post("valuesranks_enrichment_submit", {
            "api_key":         api_key,
            "identifiers":     identifiers_str,
            "species":         species,
            "ge_fdr":          fdr_threshold,
            "ge_enrichment_rank_direction": rank_direction,
            "caller_identity": caller_id,
        })
        r = result[0]
        if r.get("status") == "error":
            msg = r.get("message", "")
            # Retry on key-activation errors
            if any(kw in msg.lower() for kw in ("key", "auth", "invalid",
                                                  "not found", "not activated")):
                elapsed = time.time() - t0
                if elapsed < max_auth_wait:
                    log.info(f"  API key not yet active ({msg}); "
                             f"retrying in 30 s ({elapsed:.0f} s elapsed)…")
                    time.sleep(30)
                    continue
            raise RuntimeError(f"Submit error: {msg}")
        return r["job_id"]


# ── Poll all jobs ────────────────────────────────────────────────────

# Delay between individual status-check requests.  STRING requires ≥1 s
# between calls.  At 1.5 s/check × 100 jobs each round fires at ~0.67 req/s,
# safely within the limit.  The previous value of 0.25 s caused 429 errors.
_STATUS_POLL_DELAY = 1.5    # seconds  (≥1 s required by STRING)

# Maximum consecutive *non-429* errors before declaring a job failed.
# 429 rate-limit errors are handled separately and do NOT count here.
# Raised from 10 → 20 to tolerate transient server slowness.
_MAX_CONSECUTIVE_400 = 20

def poll_jobs(api_key, pending):
    """Poll ``{pheno: job_id}`` until every job succeeds or fails.

    Returns ``{pheno: status_dict}`` for **successful** jobs.
    """
    completed = {}
    failed = []
    remaining = dict(pending)
    # Track consecutive 400/error counts per job
    error_counts = {p: 0 for p in remaining}
    t0 = time.time()

    while remaining:
        elapsed = time.time() - t0
        if elapsed > max_wait:
            for p in remaining:
                log.error(f"  {p}: timed out after {max_wait}s")
                failed.append(p)
            break

        round_start = time.time()
        still_pending = {}
        n_new_done = 0
        n_new_fail = 0
        status_counts = {}  # status_string → count
        got_rate_limit = False  # set True if any 429 seen this round

        for pheno, job_id in remaining.items():
            try:
                # Don't retry 4xx on status checks — let poll_jobs
                # handle them via the consecutive-error counter
                status = api_get("valuesranks_enrichment_status", {
                    "api_key": api_key,
                    "job_id": job_id,
                }, politeness_delay=_STATUS_POLL_DELAY,
                   retry_on_client_error=False,
                   retries=1)
                r = status[0]
                st  = r.get("status", "unknown")
                msg = r.get("message", "")

                # Reset error counter on any successful response
                error_counts[pheno] = 0

                # Terminal statuses — STRING API uses several spellings
                _TERMINAL_FAIL = {"failed", "failure", "error"}
                if st == "success":
                    completed[pheno] = r
                    n_new_done += 1
                elif st in _TERMINAL_FAIL:
                    log.error(f"  {pheno} (job {job_id}): {st} — {msg}")
                    failed.append(pheno)
                    n_new_fail += 1
                else:
                    still_pending[pheno] = job_id
                    # Aggregate status for summary logging
                    status_counts[st] = status_counts.get(st, 0) + 1
            except Exception as exc:
                is_429 = (
                    isinstance(exc, requests.exceptions.HTTPError)
                    and exc.response is not None
                    and exc.response.status_code == 429
                )
                if is_429:
                    # Rate-limit from STRING — do NOT penalise the job.
                    # One 429 contaminates the whole round; record the flag
                    # and continue so remaining jobs still get checked.
                    got_rate_limit = True
                    still_pending[pheno] = job_id
                    status_counts["rate_limited"] = \
                        status_counts.get("rate_limited", 0) + 1
                else:
                    error_counts[pheno] = error_counts.get(pheno, 0) + 1
                    n_err = error_counts[pheno]
                    if n_err >= _MAX_CONSECUTIVE_400:
                        log.error(
                            f"  {pheno} (job {job_id}): dropping after "
                            f"{n_err} consecutive poll errors: {exc}"
                        )
                        failed.append(pheno)
                        n_new_fail += 1
                    else:
                        # Keep it pending; log only occasionally to reduce noise
                        if n_err <= 2 or n_err % 5 == 0:
                            log.warning(
                                f"  {pheno}: poll error #{n_err}: {exc}"
                            )
                        still_pending[pheno] = job_id
                        status_counts["poll_error"] = \
                            status_counts.get("poll_error", 0) + 1

        # If the round was rate-limited, back off before sleeping poll_interval.
        if got_rate_limit:
            log.warning(
                f"  Rate limited (429) during poll round — "
                f"backing off 60s before next round"
            )
            time.sleep(60)

        remaining = still_pending
        round_elapsed = time.time() - round_start
        elapsed = time.time() - t0

        if remaining:
            # Log a compact summary instead of per-job lines
            status_summary = ", ".join(
                f"{cnt} {st}" for st, cnt in sorted(status_counts.items())
            )
            log.info(
                f"  Poll round: {len(remaining)} pending "
                f"[{status_summary}], "
                f"{len(completed)} done, {len(failed)} failed "
                f"(+{n_new_done} done, +{n_new_fail} fail this round) — "
                f"round took {round_elapsed:.0f}s, "
                f"total {elapsed:.0f}s — sleeping {poll_interval}s"
            )
            time.sleep(poll_interval)

    log.info(f"Polling finished: {len(completed)} succeeded, {len(failed)} failed")
    if failed:
        log.info(f"  Failed phenotypes: {', '.join(sorted(failed))}")
    return completed


# ── Download results ─────────────────────────────────────────────────

def download_results(url):
    """Download the enrichment TSV from a STRING download_url."""
    time.sleep(1)
    resp = requests.get(url, timeout=300)
    resp.raise_for_status()
    txt = resp.text.strip()
    if not txt:
        return pd.DataFrame()
    return pd.read_csv(io.StringIO(txt), sep="\t")


# ── Fetch PPI network (edge list) ────────────────────────────────────

def fetch_ppi_network(protein_ids, score_threshold=400, batch_size=2000):
    """Retrieve the STRING PPI network among *protein_ids*.

    Calls ``/api/tsv/network`` and returns a DataFrame with columns:
      from, to, preferredName_A, preferredName_B, score, nscore,
      fscore, pscore, ascore, escore, dscore, tscore.

    This edge list loads directly into igraph::

        import igraph as ig
        g = ig.Graph.DataFrame(df[["from", "to"]], directed=False)
        # or in R:  g <- igraph::graph_from_data_frame(df, directed=FALSE)
    """
    all_edges = []
    ids = list(set(protein_ids))

    for i in range(0, len(ids), batch_size):
        batch = ids[i : i + batch_size]
        try:
            txt = api_post(
                "network",
                {
                    "identifiers": "\r".join(batch),
                    "species": species,
                    "required_score": score_threshold,
                    "caller_identity": caller_id,
                },
                fmt="tsv",
            )
            if txt.strip():
                df = pd.read_csv(io.StringIO(txt), sep="\t")
                all_edges.append(df)
        except Exception as exc:
            log.warning(f"  fetch_ppi_network batch {i}: {exc}")

    if not all_edges:
        return pd.DataFrame()
    combined = pd.concat(all_edges, ignore_index=True)

    # Standardise column names for igraph
    rename = {}
    if "preferredName_A" in combined.columns:
        rename["preferredName_A"] = "from"
    if "preferredName_B" in combined.columns:
        rename["preferredName_B"] = "to"
    combined.rename(columns=rename, inplace=True)
    return combined


# ── Build term–gene bipartite edge list ──────────────────────────────

def build_bipartite_edges(enrich_df):
    """Explode enrichment results into a term↔gene edge list.

    Returns a DataFrame with columns:
      term_id, term_name, source, gene, phenotype, fdr

    Load into igraph as a bipartite graph::

        import igraph as ig
        nodes = list(set(df['term_id']) | set(df['gene']))
        g = ig.Graph.DataFrame(df[['term_id','gene']], directed=False)
        g.vs['type'] = [v['name'] in set(df['term_id']) for v in g.vs]
    """
    rows = []
    gene_col = (
        "intersection_genes" if "intersection_genes" in enrich_df.columns
        else "query_gene_names" if "query_gene_names" in enrich_df.columns
        else None
    )
    if gene_col is None:
        return pd.DataFrame()

    for _, row in enrich_df.iterrows():
        genes_str = str(row.get(gene_col, ""))
        if not genes_str or genes_str == "nan":
            continue
        genes = [g.strip() for g in genes_str.replace(";", ",").split(",")
                 if g.strip()]
        for gene in genes:
            rows.append({
                "term_id":   row.get("term_id", ""),
                "term_name": row.get("term_name", ""),
                "source":    row.get("source", ""),
                "gene":      gene,
                "phenotype": row.get("phenotype", ""),
                "fdr":       row.get("fdr", ""),
            })
    return pd.DataFrame(rows)


# ── Column renaming for pipeline compatibility ───────────────────────
# The STRING API returns camelCase column names.  Map them to
# the snake_case / short names expected by downstream code and
# by build_bipartite_edges().
RENAME = {
    # --- current STRING API camelCase names ---
    "category":             "source",
    "termID":               "term_id",
    "termDescription":      "term_name",
    "genesMapped":          "intersection_size",
    "genesInSet":           "term_size",
    "enrichmentScore":      "enrichment_score",
    "direction":            "direction",
    "falseDiscoveryRate":   "fdr",
    "method":               "test_method",
    "proteinIDs":           "protein_ids",
    "proteinLabels":        "intersection_genes",
    "proteinInputLabels":   "query_gene_names",
    "proteinInputValues":   "input_values",
    "proteinRanks":         "protein_ranks",
    # --- legacy spaced names (older STRING versions) ---
    "Category":             "source",
    "Term id":              "term_id",
    "Term description":     "term_name",
    "Genes mapped":         "intersection_size",
    "Genes in set":         "term_size",
    "Enrichment score":     "enrichment_score",
    "Direction":            "direction",
    "Count in set/pathway": "count_in_set",
    "False Discovery Rate": "fdr",
    "Method":               "test_method",
}


# ══════════════════════════════════════════════════════════════════════
# Main execution
# ══════════════════════════════════════════════════════════════════════

t_start = time.time()

# 1. API key
api_key = get_api_key()
log.info(f"API key: {api_key[:8]}…")

# 2. Discover phenotype files
gene_set_dir = Path(str(INPUT.dir_gene_set))
files = {f.stem: f for f in sorted(gene_set_dir.glob("*.txt"))}
log.info(f"Found {len(files)} phenotype file(s) in {gene_set_dir}")

# 3. (Optional) pre-map universe to STRING IDs for faster processing
id_map = {}
if premap_ids:
    universe = [l.strip() for l in open(str(INPUT.gene_universe)) if l.strip()]
    log.info(f"Pre-mapping {len(universe)} universe genes → STRING IDs …")
    id_map = map_identifiers(universe)
    log.info(f"  Mapped {len(id_map)} / {len(universe)} genes")

# 4. Read files & submit jobs
pending = {}   # pheno → job_id
skipped = []
submit_failed = {}  # pheno → ident_str  (for retry passes)
_n_submitted_this_batch = 0  # counts successfully submitted jobs in current batch

for i, (pheno, fpath) in enumerate(files.items(), 1):
    log.info(f"[{i}/{len(files)}] {pheno}: reading {fpath}")
    ident_str, n, ranked = read_gene_file(fpath)

    if not ranked:
        log.info(f"[{i}/{len(files)}] {pheno}: "
                 "not a 2-column (gene<TAB>value) file — skipping")
        skipped.append(pheno)
        continue

    if n < min_genes:
        log.info(f"[{i}/{len(files)}] {pheno}: "
                 f"{n} genes < min_genes={min_genes} — skipping")
        skipped.append(pheno)
        continue

    # Optional: replace gene names with STRING IDs
    if premap_ids and id_map:
        new_lines = []
        for line in ident_str.split("\n"):
            gene, val = line.split("\t", 1)
            new_lines.append(f"{id_map.get(gene, gene)}\t{val}")
        ident_str = "\n".join(new_lines)

    log.info(f"[{i}/{len(files)}] {pheno}: {n} genes — submitting …")
    try:
        jid = submit_job(ident_str, api_key)
        pending[pheno] = jid
        _n_submitted_this_batch += 1
        log.info(f"[{i}/{len(files)}] {pheno}: job_id = {jid}")
        # Batch pause: after every submit_batch_size successful submissions,
        # sleep to avoid flooding STRING's queue with too many large jobs at
        # once (server-side 'globalenrichment.py error code 1' failures).
        if submit_batch_size > 0 and _n_submitted_this_batch % submit_batch_size == 0:
            log.info(
                f"  Batch of {submit_batch_size} submitted — "
                f"pausing {submit_batch_pause}s before next batch …"
            )
            time.sleep(submit_batch_pause)
    except Exception as exc:
        log.error(f"[{i}/{len(files)}] {pheno}: submit ERROR — {exc}")
        submit_failed[pheno] = ident_str

# 4b. Retry failed submissions (502s are transient server errors)
#      Wait with increasing backoff; probe server health before each pass
#      to avoid wasting retries while server is still down.
_SUBMIT_RETRY_WAITS = [60, 120, 180, 300, 600]  # seconds between passes

def _server_is_healthy(max_probe_wait=300):
    """Probe STRING API until it responds or *max_probe_wait* elapses."""
    t0 = time.time()
    probe_interval = 30
    while time.time() - t0 < max_probe_wait:
        try:
            r = requests.get(
                f"{STRING_API}/json/version", timeout=30
            )
            if r.status_code < 500:
                return True
        except Exception:
            pass
        log.info(f"  Server not healthy yet — retrying in {probe_interval}s …")
        time.sleep(probe_interval)
    return False

for retry_pass, wait in enumerate(_SUBMIT_RETRY_WAITS, 1):
    if not submit_failed:
        break
    log.info(
        f"Submit retry pass {retry_pass}/{len(_SUBMIT_RETRY_WAITS)}: "
        f"{len(submit_failed)} failed submission(s) — "
        f"waiting {wait}s before retrying …"
    )
    time.sleep(wait)

    # Probe server health before burning through retries
    if not _server_is_healthy():
        log.warning(
            f"  Server still unhealthy after probe — "
            f"will attempt submissions anyway"
        )

    still_failed = {}
    recovered = 0
    for pheno, ident_str in submit_failed.items():
        try:
            jid = submit_job(ident_str, api_key)
            pending[pheno] = jid
            recovered += 1
            log.info(f"  {pheno}: retry OK — job_id = {jid}")
        except Exception as exc:
            log.warning(f"  {pheno}: retry pass {retry_pass} failed — {exc}")
            still_failed[pheno] = ident_str
    log.info(
        f"  Retry pass {retry_pass} result: "
        f"{recovered} recovered, {len(still_failed)} still failing"
    )
    submit_failed = still_failed

# Any phenotypes that still failed after all retry passes
if submit_failed:
    log.error(
        f"{len(submit_failed)} phenotypes permanently failed after "
        f"{len(_SUBMIT_RETRY_WAITS)} retry passes: "
        f"{', '.join(sorted(submit_failed.keys()))}"
    )
    skipped.extend(submit_failed.keys())

log.info(f"Submitted {len(pending)} job(s), skipped {len(skipped)}")

# 5. Poll all jobs
if pending:
    completed = poll_jobs(api_key, pending)
else:
    completed = {}

# 6. Download & combine results
all_dfs = []

for pheno, status in completed.items():
    dl_url = status.get("download_url")
    if not dl_url:
        log.warning(f"  {pheno}: no download_url in status — skipping")
        continue

    log.info(f"  {pheno}: downloading results …")
    try:
        df = download_results(dl_url)
    except Exception as exc:
        log.error(f"  {pheno}: download error: {exc}")
        continue

    if df.empty:
        log.info(f"  {pheno}: empty result")
        continue

    # Rename columns for pipeline compatibility
    df.rename(
        columns={k: v for k, v in RENAME.items() if k in df.columns},
        inplace=True,
    )
    df["phenotype"]   = pheno
    df["job_id"]      = pending[pheno]
    df["page_url"]    = status.get("page_url", "")
    df["significant"] = True   # STRING returns only significant terms by default

    fdr_col = "fdr" if "fdr" in df.columns else "False Discovery Rate"
    if fdr_col in df.columns:
        n_sig = (pd.to_numeric(df[fdr_col], errors="coerce") < fdr_threshold).sum()
    else:
        n_sig = len(df)
    log.info(f"  {pheno}: {len(df)} terms, {n_sig} with FDR < {fdr_threshold}")
    all_dfs.append(df)

# 7. Write enrichment output
out_path = str(OUTPUT.enrich)

if all_dfs:
    combined = pd.concat(all_dfs, ignore_index=True)
    log.info(
        f"Writing {len(combined)} rows from "
        f"{len(all_dfs)}/{len(files)} phenotypes → {out_path}"
    )
    combined.to_csv(out_path, sep="\t", index=False)
else:
    combined = pd.DataFrame()
    log.warning("No enrichment results — writing empty output")
    combined.to_csv(out_path, sep="\t", index=False)

# 8. Build & write PPI network edge list
ppi_path = str(OUTPUT.ppi_network)
if not combined.empty:
    # Collect all enriched protein labels across phenotypes
    gene_col = (
        "intersection_genes" if "intersection_genes" in combined.columns
        else "query_gene_names" if "query_gene_names" in combined.columns
        else None
    )
    all_proteins = set()
    if gene_col:
        for val in combined[gene_col].dropna():
            for g in str(val).replace(";", ",").split(","):
                g = g.strip()
                if g:
                    all_proteins.add(g)

    score_threshold = int(getattr(PARAMS, "score_threshold", 400))
    if all_proteins:
        log.info(f"Fetching PPI network for {len(all_proteins)} enriched proteins …")
        ppi_df = fetch_ppi_network(list(all_proteins),
                                   score_threshold=score_threshold)
        if not ppi_df.empty:
            log.info(f"  PPI network: {len(ppi_df)} edges")
            ppi_df.to_csv(ppi_path, sep="\t", index=False)
        else:
            log.info("  PPI network: no edges returned")
            pd.DataFrame(columns=["from", "to", "score"]).to_csv(
                ppi_path, sep="\t", index=False)
    else:
        log.info("  No protein labels found — writing empty PPI edge list")
        pd.DataFrame(columns=["from", "to", "score"]).to_csv(
            ppi_path, sep="\t", index=False)
else:
    pd.DataFrame(columns=["from", "to", "score"]).to_csv(
        ppi_path, sep="\t", index=False)
log.info(f"PPI network edge list: {ppi_path}")

# 9. Build & write term–gene bipartite edge list
bip_path = str(OUTPUT.bipartite)
if not combined.empty:
    bip_df = build_bipartite_edges(combined)
    if not bip_df.empty:
        log.info(f"Term–gene bipartite graph: {len(bip_df)} edges, "
                 f"{bip_df['term_id'].nunique()} terms, "
                 f"{bip_df['gene'].nunique()} genes")
        bip_df.to_csv(bip_path, sep="\t", index=False)
    else:
        log.info("  No bipartite edges — writing empty file")
        pd.DataFrame(columns=["term_id", "term_name", "source", "gene",
                               "phenotype", "fdr"]).to_csv(
            bip_path, sep="\t", index=False)
else:
    pd.DataFrame(columns=["term_id", "term_name", "source", "gene",
                           "phenotype", "fdr"]).to_csv(
        bip_path, sep="\t", index=False)
log.info(f"Term–gene bipartite edge list: {bip_path}")

elapsed = time.time() - t_start
log.info(f"Done ({elapsed:.1f}s)")
