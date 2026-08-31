#!/usr/bin/env bash
# Build the anonymised artifact §5 promises, ready to upload.
#
# WHAT THIS DOES NOT DO: upload it. Creating an anonymous.4open.science mirror
# needs a GitHub OAuth login, and publishing the code is an author's decision,
# not a script's. This produces the directory; you publish it and paste the URL
# into \artifacturl in main.tex.
#
# WHY A SCRUB IS NEEDED AT ALL. The repository de-anonymises the submission three
# ways, and only the third is obvious:
#
#   1. Git history. Six committer identities, two of them on an institutional
#      domain. An affiliation in an email address is the strongest de-anonymiser
#      there is, and it survives every scrub that only edits files. This script
#      copies the working tree WITHOUT .git for exactly that reason.
#   2. .github/CODEOWNERS names a GitHub handle seven times.
#   3. CONTRIBUTING.md carries the clone URL.
#
# It also enforces the release policy §5 states: derived summaries and code, no
# cell-level matrices. data/ is gitignored, so it is excluded by construction —
# but the check is explicit rather than assumed, because "it is gitignored" is
# the kind of assumption this paper is about.
set -euo pipefail
cd "$(dirname "$0")/../.."          # repo root
OUT="${1:-artifact}"

if [ -e "$OUT" ]; then
  echo "refusing to overwrite existing $OUT — remove it or pass another path" >&2
  exit 1
fi

echo "building anonymised artifact in $OUT/"

# Tracked files only, so anything gitignored (data/, caches) cannot leak. And no
# .git: history is the thing that cannot be scrubbed in place.
mkdir -p "$OUT"
git ls-files -z | while IFS= read -r -d '' f; do
  case "$f" in
    .github/CODEOWNERS) continue ;;   # GitHub handles; rewritten below
    paper/*)            continue ;;   # the submission itself does not ship
  esac
  mkdir -p "$OUT/$(dirname "$f")"
  cp "$f" "$OUT/$f"
done

# Replace, rather than drop, so a reader can see the review policy existed.
mkdir -p "$OUT/.github"
cat > "$OUT/.github/CODEOWNERS" <<'EOF'
# Redacted for double-blind review. The original assigns each frozen file
# (src/schema.py, config/panel.yaml, config/labeling_axes.yaml,
# tests/test_freeze.py, CLAUDE.md) a code owner, so CLAUDE.md invariant 3 —
# two approvals and a written reason — is enforced by branch protection rather
# than by convention.
EOF

# Scrub the identifying strings that survive in file content.
python3 - "$OUT" <<'PY'
import pathlib, re, sys
root = pathlib.Path(sys.argv[1])
SUBS = [
    (re.compile(r"https://github\.com/[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+(\.git)?"),
     "https://anonymous.4open.science/r/ANONYMISED"),
    (re.compile(r"@diegomardian\b"), "@owner"),
    (re.compile(r"\bdiegomardian\b"), "anonymised"),
    (re.compile(r"\bMeowsers4\b|\bbodemeowsers\b|\bbodebosell\b", re.I), "anonymised"),
    (re.compile(r"/Users/[A-Za-z0-9_.-]+/"), "/path/to/"),
]
TEXT = {".md", ".py", ".yml", ".yaml", ".toml", ".cfg", ".txt", ".sh", ".json"}
changed = 0
for p in root.rglob("*"):
    if not p.is_file() or p.suffix.lower() not in TEXT:
        continue
    try:
        s = p.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        continue
    new = s
    for pat, rep in SUBS:
        new = pat.sub(rep, new)
    if new != s:
        p.write_text(new, encoding="utf-8")
        changed += 1
print(f"  scrubbed {changed} files")
PY

echo
echo "verifying:"
fail=0
if [ -d "$OUT/.git" ]; then echo "  FAIL  .git present"; fail=1
else echo "  ok    no .git — history cannot leak"; fi

if find "$OUT" -type f \( -name '*.h5ad' -o -name '*.h5' -o -name '*.loom' \
     -o -name '*.mtx' -o -name '*.mtx.gz' -o -name '*.rds' -o -name '*.RData' \) \
     | grep -q .; then
  echo "  FAIL  cell-level matrices present — §5 promises derived summaries only"; fail=1
else
  echo "  ok    no cell-level matrices"
fi

if hits=$(grep -rniE 'diegomardian|bodebosell|meowsers|github\.com/[A-Za-z0-9_-]+|@bu\.edu|/Users/' \
          "$OUT" 2>/dev/null); then
  echo "  FAIL  identifying strings survive the scrub:"; echo "$hits" | sed 's/^/        /' | head -20
  fail=1
else
  echo "  ok    no identifying strings in any shipped file"
fi

echo "  ..    $(find "$OUT" -type f | wc -l | tr -d ' ') files, $(du -sh "$OUT" | cut -f1)"
echo
if [ "$fail" -eq 0 ]; then
  cat <<EOF
PASS. Next, and only you can do these:

  1. Publish $OUT/ as an anonymised mirror (anonymous.4open.science, or a
     fresh repo with no history under a throwaway account).
  2. Paste the URL into \\artifacturl in paper/wmhs/main.tex.
  3. Rebuild, then run paper/wmhs/check_anonymity.sh — it must exit 0.
EOF
else
  echo "NOT SAFE TO PUBLISH — fix the failures above."
fi
exit "$fail"
