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
| Full paper | `main.tex` | **8 pages, measured** | 9 |
| Extended abstract | `main_short.tex` | **4 pages, measured** | 4 |

**Measured 2026-09-11 against real NeurIPS geometry, and both pass
`./build.sh`.** `neurips_2026.sty` is still not published, so the build was run
against the official `neurips_2024.sty` (from
`media.neurips.cc/Conferences/NeurIPS2024/Styles.zip`) behind a five-line shim
absorbing `[dblblindworkshop]` and `\workshoptitle`. NeurIPS page geometry has
been stable across these years, so this measures pages accurately even though it
is **not** the file to submit with.

**The previous "9 pages, measured" was measured on a broken document, and that
is this paper's own thesis arriving in its own page check.** Fifteen table rows
across `sections/calibration.tex` and `sections/appendix.tex` ended in a single
`\` instead of `\\` — a control space, not a row break — so the three-grid table
and the residual matrix collapsed. `pdflatex` emitted 17 errors and then hung;
TeX dropped the table bodies, and the page count taken off that build came in
about a page short. A literal TAB had also eaten the `\t` of `\times` in
`sections/conclusion.tex`, rendering "0.82imes". All fixed. **A page count is
only a measurement if `errors` is 0 — read that line before the page line.**

With the errors fixed the full build ran to ten pages of main text, so §4, the
single-cell competitor benchmark, moved to Appendix A in **both** builds. It is
the most instance-specific section in the paper and the one neither build's
argument depends on; the conclusion now carries its result and a pointer in both
builds. The full build is 8 pages with room to spare.

The four-page build was also over, by about eleven lines, and had been before
any of this. It is in now by: a tightened abstract, the related-work paragraph
pointing at the appendix survey both builds already carry, Figure 1 at 0.88
width, and four restatements gated to the full build. Every cut is prose the
short build's own appendix carries verbatim or that the same build states
elsewhere. No number moved.

Neither the shim nor `neurips_2024.sty` is committed. **Download the real 2026
workshop style when it is published and re-run `./build.sh` before submitting**
--- this measurement is a floor on confidence, not a substitute.

Reproduce:

```
curl -L -o /tmp/s.zip https://media.neurips.cc/Conferences/NeurIPS2024/Styles.zip
unzip -j /tmp/s.zip 'Styles/neurips_2024.sty' -d paper/wmhs/
printf '%s\n' '\NeedsTeXFormat{LaTeX2e}' '\ProvidesPackage{neurips_2026}' \
  '\DeclareOption{dblblindworkshop}{}' \
  '\DeclareOption*{\PassOptionsToPackage{\CurrentOption}{neurips_2024}}' \
  '\ProcessOptions\relax' '\RequirePackage{neurips_2024}' \
  '\newcommand{\workshoptitle}[1]{}' > paper/wmhs/neurips_2026.sty
cd paper/wmhs && ./build.sh
```

Limits verified against the CFP's own text on 2026-08-31: *"Full Papers: at most
9 pages of main text. Extended Abstracts: at most 4 pages of main text.
References and appendices do not count toward the page limit, but the main text
must be self-contained."* The extended-abstract cap is **4**, not 6.

Both `\input` the same `sections/` files and differ only in the `\iffull` flag
set at the top. **The short build is a strict subset of the same prose, not a
rewrite** — so a number cannot say one thing in one version and something else
in the other, which given what this paper argues is the one way it must not be
wrong. Whatever the short build drops from the main text, its appendix carries
verbatim. The one place the two files hold separate prose is the abstract, which
lives in `main.tex` and `main_short.tex` rather than in `sections/`; the short
one states fewer results, never different ones.

The short build drops from the main text: the sound-and-not-complete analysis
with the residual matrix and the information ratio, the censoring experiment,
the OLS-splits-nothing argument, the second cohort, the three grids, and four
restatements. It keeps the whole of §2 — the general form, the
trial-simulator instantiation and the one-line check — Figure 1, and the
mandatory responsible-use statement.

**If you are choosing: submit `main.tex`.** The material only its main text
carries — the residual matrix, the information ratio, and the censoring
result — is what answers the "this is just the inverse crime" objection and
is the most venue-specific evidence in the paper. The short build is a good
fallback and its appendix carries all of it, but a reviewer who forms a novelty
judgment on page 2 never gets there.

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
python paper/wmhs/make_fig3.py     # Figure 1 — the generator statistic
python paper/wmhs/make_fig4.py     # Figure 2 — the trial recovery curve
python paper/wmhs/make_fig1.py     # Figure 3 — the calibration
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
- [x] ~~Re-check the page limits.~~ Done 2026-09-11 against the real geometry:
      full 8 of 9, short 4 of 4, `./build.sh` PASS with `errors 0`. Re-run it
      after any edit — both builds now sit inside their limits, the short one
      exactly at its.
- [ ] **Read the `errors` line before the page line.** A build with LaTeX errors
      drops content and under-reports its own length; that is how the previous
      "9 pages, measured" happened. `./build.sh` fails on `errors > 0`, on
      undefined references, and now on a missing `pypdf` rather than skipping
      the page check.
- [ ] Rebuild both figures from freshly written tables and confirm the captions
      still match what the tables say. `pytest tests/test_paper_numbers.py`
      checks the prose against the tables it quotes; it is not a substitute for
      looking at the figures.
- [ ] `grep -c 'undefined' main.log` returns 0 after a full four-pass build.
- [x] ~~Verify the `petukhov2022` (cacoa) venue.~~ Done 2026-08-31: cacoa is
      still a preprint, no journal version exists, so the bioRxiv entry is
      correct. All eight added entries were checked against the publishers'
      own records — see the header of `refs.bib`.
- [ ] Decide which to submit. Both build and both fit. **`main.tex` is the
      recommendation** and the reasons are under "Two builds, one source";
      `main_short.tex` is ready if you want it.

## Layout

```
main.tex                 full paper — preamble, abstract, \input list, \fulltrue
main_short.tex           extended abstract — the same, \fullfalse
sections/setup.tex       §1  the decomposition and the rule
sections/blind.tex       §2  the recovery curve, the general form, the
                             trial-simulator instantiation, the one-line check
sections/calibration.tex §3  what a correct calibration returns (Figure 3)
sections/withdrawn.tex       three more checks that could not fire — appendix
sections/conclusion.tex  §4  the closer, and the pointer to Appendix A
sections/responsible.tex §5  the mandatory responsible-use statement
sections/bench.tex           Appendix A — what abstention buys, BOTH builds
sections/benchtable.tex      Table 1, appendix in both builds
sections/gridcaveats.tex     appendix in both builds
sections/modaloutcome.tex    ditto
sections/relatedwork.tex     ditto
sections/trialtable.tex      Table 2, appendix in both builds
sections/appendix.tex    the appendices, and Table 3 (the residual matrix)
sections/llm.tex         final appendix, LLM use
refs.bib                bibliography
make_fig1.py            Figure 3, the calibration
make_fig3.py            Figure 1, the generator statistic
make_fig4.py            Figure 2, the trial recovery curve
_tables.py              resolves a result table by name
build.sh                builds both and enforces both page limits
make_overleaf.sh        bundles overleaf.zip (tex + bib + figures only)
check_anonymity.sh      double-blind guard; run before the final build
make_artifact.sh        builds an anonymised code release, if one is ever wanted
PATCHES.md              working record of corrections, with what is still open
```
