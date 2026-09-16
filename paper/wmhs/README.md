# WMHS @ NeurIPS 2026 — submission source

## Current full-paper revision — 2026-09-13

`main.tex` reads `sections/full/`, including the revised abstract and condensed
appendix. The extended abstract reads `sections/`; its one stale historical-grid
qualification was updated without adding the new full-paper results. The full paper uses an editorial
variant of Figure 3 labelled "candidate", regenerated from the existing table
with `python paper/wmhs/make_fig1_full.py` from the repo root; the original
figure used by the shorter version is unchanged.

`results_manifest.json` pins every result table read by the paper checks and
figure scripts. Re-running an experiment no longer changes the submission by
directory ordering; update the manifest explicitly when adopting a new result.
CI requires every pinned file, follows the recursive `main.tex` input graph
through `sections/full/`, and compiles both manuscript variants.

A follow-up pass adds the saved Bernoulli-randomisation and RMST controls, a
simulation-settings table, and the final controlled-grid, inner-budget and
seed-preserving patient-influence results. Untraceable historical rates and
repetitive diagnostic history were removed. The completed follow-up plan is in
[the additional-runs plan](../../docs/wmhs_next_runs_plan.md).
Fresh-environment commands for the independent lifelines comparison and the
residual/performance clean control are in [REPRODUCE.md](REPRODUCE.md).

The final full-paper review promotes the saved clean-control comparison to a
main-text table, adds controlled-grid and fixed-pair sensitivity tables, and
reports the independent Cox implementation check. The main text now specifies
the calibration substrate and interval sampling units, distinguishes low cell
counts from an absent population, and qualifies the untested learned-generator
implications. Per-count candidates and unobserved lower crossing boundaries
make the limits of stable binned cutpoints explicit. No new experiments were
run for this editorial expansion; the appendix and shorter paper are unchanged.

Both manuscripts are compiled and page-limited in CI. The full build's main
text fills 9 pages (limit 9); references start on page 10 and the condensed
appendix on page 12, for 15 pages total. The extended abstract is unchanged
at 4 pages of main text (limit 4). The full build has
zero LaTeX errors, undefined references/citations, or overfull boxes. The
anonymity check passes. Paper-number tests now cover claims in both source
trees, including the newly adopted controlled and sensitivity results.

The tracked `neurips_2026.sty` is the local compatibility shim over the
2024 style described below. These page counts apply to that installed style.
The official submission style still needs to be supplied before submission.

Build the full paper alone from this directory:

```sh
pdflatex -halt-on-error main.tex
bibtex main
pdflatex -halt-on-error main.tex
pdflatex -halt-on-error main.tex
./check_anonymity.sh
```

`build.sh` still builds both versions. `make_overleaf.sh` now includes the
nested full-paper sources. The notes below describe the earlier shared-source
build and are retained as historical context; their page counts, section map
and statements that both versions use identical prose are superseded above.

---


**Venue:** [World Models for High-Stakes Health: Reliable Clinical Trial
Simulation and Intervention-Aware Reasoning](https://wmhs-neurips.github.io/WMHS/),
NeurIPS 2026, Atlanta.
**Deadline:** 15 September 2026, AoE. **Non-archival, double-blind.**
**Limits:** 9 pages of main text (full paper) or 4 (extended abstract).
References and appendices do not count. Verified against the CFP's own text
on 2026-08-31: *"Full Papers: at most 9 pages of main text. Extended
Abstracts: at most 4 pages of main text. References and appendices do not
count toward the page limit, but the main text must be self-contained."*
A reviewer could not locate this limit and raised over-length as a risk, so
the source is quoted here rather than paraphrased.
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
single-cell competitor benchmark, moved into the appendix in **both** builds. It is
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

The current builds have distinct source trees: `main.tex` reads
`sections/full/`, while `main_short.tex` reads `sections/`. Shared numerical
claims are checked against both source files and their pinned result tables.

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

The tracked `neurips_2026.sty` is a compatibility shim over the tracked official
2024 style. Replace it with the official workshop style when that file is
published, rerun the build gates, and do not modify the official style.

## Figures — regenerate, never transcribe

```
python paper/wmhs/make_fig3.py     # Figure 1 — the generator statistic
python paper/wmhs/make_fig4.py     # Figure 2 — the trial recovery curve
python paper/wmhs/make_fig1.py     # Figure 3 — the calibration
```

The scripts read paths pinned in `results_manifest.json` and **print the path
they used**, so every number on an axis traces to a versioned table with a
commit hash and a fixed seed. A later result cannot silently change a figure.

The tables come from:

```
python -m src.harness.calibration_gap                     # 50 replicates, 13 seeds
python -m src.harness.calibration_gap --replicates 500 --seeds 8
```

That module is the sweep behind Figure 1. It did not exist when the figure was
first drawn — the run came from an uncommitted script, which is item 3 of the
implementation-defects appendix, and re-deriving it is why the module is here.

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

## Appendix order

Reordered 2026-09-11 to follow the main text rather than the order things were
written in. §2's support first, then §3's, then the benchmark as one run of
three consecutive sections instead of three scattered ones, then context and
record:

```
A  The curve itself                    G  What abstention buys, and which part
B  Five estimators of the same effect  H  The competitor benchmark
C  The residual matrix                 I  What abstention buys, per world
D  What the grid was reporting         J  Related work
E  Abstention is the modal outcome     K  Three implementation defects
F  Three guards that could not fire    L  Extended limitations
                                       M  Use of large language models
```

F and K used to be called "Three guards that could not fire" and "The withdrawn
guards" — two titles for two different things. F is the analysis guards, K is
the implementation defects.

## Layout

```
main.tex                 full paper — preamble, abstract, \input list, \fulltrue
main_short.tex           extended abstract — the same, \fullfalse
sections/setup.tex       §1  the decomposition and the rule
sections/blind.tex       §2  the recovery curve, the general form, the
                             trial-simulator instantiation, the one-line check
sections/calibration.tex §3  what a correct calibration returns (Figure 3)
sections/withdrawn.tex       three more checks that could not fire — appendix
sections/conclusion.tex  §4  the closer, and the pointer to the benchmark
sections/responsible.tex §5  the mandatory responsible-use statement
sections/bench.tex           appendix in BOTH builds — what abstention buys
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
