# WMHS @ NeurIPS 2026 — submission source

**Venue:** [World Models for High-Stakes Health: Reliable Clinical Trial
Simulation and Intervention-Aware Reasoning](https://wmhs-neurips.github.io/WMHS/),
NeurIPS 2026, Atlanta.
**Deadline:** 15 September 2026, AoE. **Non-archival, double-blind.**
**Limits:** 9 pages of main text (full paper) or 4 (extended abstract).
References and appendices do not count. Both versions are built and both fit.

A **responsible-use statement covering limitations and impact is mandatory** —
a submission without one is desk-rejected. It is `sections/responsible.tex`,
rendered as §5. Do not cut it for space.

## Two builds, one source

| | file | main text | limit |
|---|---|---|---|
| Full paper | `main.tex` | 7 pages | 9 |
| Extended abstract | `main_short.tex` | 4 pages | 4 |

Limits verified against the CFP's own text on 2026-08-31: *"Full Papers: at most
9 pages of main text. Extended Abstracts: at most 4 pages of main text.
References and appendices do not count toward the page limit, but the main text
must be self-contained."* The extended-abstract cap is **4**, not 6.

Both `\input` the same `sections/` files and differ only in the `\iffull` flag
set at the top. **The short build is a strict subset of the same prose, not a
rewrite** — so a number cannot say one thing in one version and something else
in the other, which given what this paper argues is the one way it must not be
wrong. Whatever the short build drops from the main text, its appendix carries
verbatim.

The short build drops: the related-work survey, Table 1 and the invented-effect
detail, the conformal aside, the grid caveats, the modal-outcome result, and the
two withdrawn guards. It keeps the whole of §2 — the general form, the
trial-simulator instantiation and the one-line check — Figure 1, and the
mandatory responsible-use statement.

```
./build.sh          # builds both and FAILS if either exceeds its page limit
./make_overleaf.sh  # bundles overleaf.zip for upload
```

Or by hand: `pdflatex main && bibtex main && pdflatex main && pdflatex main`.

`neurips_2026.sty` is **not vendored here.** Download the official workshop
style from the NeurIPS site and drop it beside `main.tex`. The document loads it
as `\usepackage[dblblindworkshop]{neurips_2026}` with `\workshoptitle{...}`,
which is the workshop track's own interface — do not modify the style file, the
style file says tweaking it risks desk rejection, and the page budget is met by
writing less.

## Figures — regenerate, never transcribe

```
python paper/wmhs/make_fig1.py     # Figure 1  — the calibration
python paper/wmhs/make_fig3.py     # Figure 2  — the recovery curve
```

Both read the newest matching table under `results/` and **print the path they
used**, so every number on an axis traces to a versioned table with a commit
hash and a fixed seed. Neither hard-codes a sha: the sweep gets re-derived, and
a hard-coded path silently goes stale rather than failing.

The tables come from:

```
python -m src.harness.calibration_gap                     # 50 replicates, 13 seeds
python -m src.harness.calibration_gap --replicates 500 --seeds 8
```

That module is the sweep behind Figure 1. It did not exist when the figure was
first drawn — the run came from an uncommitted script, which is Appendix A
item 3, and re-deriving the number is why the module is here.

## Before submitting

- [ ] **Run `./check_anonymity.sh` after the final build.** It greps the
      sources *and* the built PDF for identifying strings and checks the PDF
      metadata, which is where an author name leaks without appearing on any
      page. Exit 0 means safe. It passes today.

      No artifact link ships: §5 states the release *policy* — derived summaries
      and code, never cell-level matrices — without pointing at a URL. If that
      ever changes, `./make_artifact.sh` builds an anonymised copy, and note
      that the repository de-anonymises three ways: **git history** carries six
      committer identities, two on an institutional domain (an affiliation in an
      email address survives every scrub that only edits files, which is why
      that script drops `.git`), `.github/CODEOWNERS` names a GitHub handle
      seven times, and `CONTRIBUTING.md` carries the clone URL.
- [ ] Rebuild both figures from freshly written tables and confirm the captions
      still match what the tables say.
- [ ] `grep -c 'undefined' main.log` returns 0 after a full four-pass build.
- [x] ~~Verify the `petukhov2022` (cacoa) venue.~~ Done 2026-08-31: cacoa is
      still a preprint, no journal version exists, so the bioRxiv entry is
      correct. All eight added entries were checked against the publishers'
      own records — see the header of `refs.bib`.
- [ ] Decide which to submit. Both build and both fit; `main.tex` is the
      recommendation, `main_short.tex` is ready if you want it.

## Layout

```
main.tex                 full paper — preamble, abstract, \input list, \fulltrue
main_short.tex           extended abstract — the same, \fullfalse
sections/setup.tex       §1  the decomposition and the rule
sections/blind.tex       §2  the recovery curve, the general form, the
                             trial-simulator instantiation, the one-line check
sections/calibration.tex §3  what a correct calibration returns (Figure 1)
sections/withdrawn.tex       two more checks that could not fire — inside §3
sections/bench.tex       §4  what abstention buys (Table 1), and the closer
sections/responsible.tex §5  the mandatory responsible-use statement
sections/benchtable.tex      Table 1 — main text (full) or appendix (short)
sections/gridcaveats.tex     ditto
sections/modaloutcome.tex    ditto
sections/appendix.tex    Appendices A-D
sections/llm.tex         Appendix E, LLM use
refs.bib                bibliography
make_fig1.py            Figure 1
make_fig3.py            Figure 2
_tables.py              resolves a result table by name
build.sh                builds both and enforces both page limits
make_overleaf.sh        bundles overleaf.zip (tex + bib + figures only)
check_anonymity.sh      double-blind guard; run before the final build
make_artifact.sh        builds an anonymised code release, if one is ever wanted
PATCHES.md              working record of corrections, with what is still open
```
