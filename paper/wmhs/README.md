# WMHS @ NeurIPS 2026 — submission source

**Venue:** [World Models for High-Stakes Health: Reliable Clinical Trial
Simulation and Intervention-Aware Reasoning](https://wmhs-neurips.github.io/WMHS/),
NeurIPS 2026, Atlanta.
**Deadline:** 15 September 2026, AoE. **Non-archival, double-blind.**
**Limits:** 9 pages of main text (full paper) or 4 (extended abstract).
References and appendices do not count. Currently ~6.5 pages of main text.

A **responsible-use statement covering limitations and impact is mandatory** —
a submission without one is desk-rejected. It is `sections/responsible.tex`,
rendered as §5. Do not cut it for space.

## Build

```
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

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

- [ ] **Fill in `\artifacturl` in `main.tex`.** §5 promises a code release and
      the submission is double-blind, so it must be an anonymising mirror
      (`anonymous.4open.science`), never `github.com/...`. Until you do,
      `./check_anonymity.sh` fails and the placeholder renders **in red** in
      the PDF, so it cannot ship unnoticed.
- [ ] **Run `./check_anonymity.sh` after the final build.** It greps the
      sources *and* the built PDF for identifying strings, checks the PDF
      metadata, and fails on the unfilled placeholder. Exit 0 means safe.
- [ ] Rebuild both figures from freshly written tables and confirm the captions
      still match what the tables say.
- [ ] `grep -c 'undefined' main.log` returns 0 after a full four-pass build.
- [x] ~~Verify the `petukhov2022` (cacoa) venue.~~ Done 2026-08-31: cacoa is
      still a preprint, no journal version exists, so the bioRxiv entry is
      correct. All eight added entries were checked against the publishers'
      own records — see the header of `refs.bib`.
- [ ] Decide full paper (9pp) vs extended abstract (4pp). The section order was
      chosen so the first four pages stand alone: §1-§3 carry the headline
      result, the general form, the venue instantiation and Figure 1. A 4-page
      version is `setup` + `blind` + `calibration` with §4's Table 1 moved to
      an appendix.

## Layout

```
main.tex                 preamble, abstract, and the \input list
sections/setup.tex       §1  the decomposition and the rule
sections/blind.tex       §2  the recovery curve, the general form, the
                             trial-simulator instantiation, the one-line check
sections/calibration.tex §3  what a correct calibration returns (Figure 1)
sections/withdrawn.tex       two more checks that could not fire — inside §3
sections/bench.tex       §4  what abstention buys (Table 1), and the closer
sections/responsible.tex §5  the mandatory responsible-use statement
sections/appendix.tex    Appendices A-D
sections/llm.tex         Appendix E, LLM use
refs.bib                bibliography
make_fig1.py            Figure 1
make_fig3.py            Figure 2
_tables.py              resolves a result table by name
check_anonymity.sh      double-blind guard; run before the final build
PATCHES.md              working record of corrections, with what is still open
```
