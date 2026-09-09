# A2: request access, then inventory metadata only

This is the next A2 action after the public-product inventory. It seeks only
to establish whether the original Chen/HTAN release contains a MANDO-derived
cell table or segmentation. It does not request a scientific result and does
not authorize downloading raw images.

**Authenticated result, 2026-09-09.** The metadata-only inventory was run with
the configured account and returned HTTP 403 for `syn23520239`:
`You lack READ access to the requested entity.` The job then stopped rather than
inventing a partial inventory or treating the inaccessible entity as absent. No
file content was downloaded and no result artifact was written. The account-only
test is therefore complete; the access request below is now necessary.

## Access request text

Submit through the Sage service desk linked from the [existing request for these
same entities](https://www.synapse.org/Synapse%3Asyn5637528/discussion/threadId%3D11828).
Replace the bracketed fields before sending:

> **Subject:** Read-only metadata access request: COLON MAP MxIF / HTA11
>
> I am [name, academic affiliation or independent-researcher status] reusing
> the published COLON MAP / Chen et al. HTAN data for a non-commercial
> reproducibility analysis. Please grant read access to `syn23564801`,
> `syn23630431`, and `syn23520239`. My Synapse username is [username].
>
> My immediate use is metadata inventory only: I will enumerate entity names,
> file types, and annotations to determine whether MANDO-derived cell-level
> intensity tables or segmentation masks exist. I will not download raw images
> or redistribute data. If a suitable processed table exists, I will use it
> only for the described same-cohort protein validation, with patient-level
> aggregation and no attempt to infer a per-cell GUCA2A result from this panel.
>
> I authenticated successfully, but a metadata-only request returns HTTP 403 for
> `syn23520239`. Please let me know whether a data-use agreement, training, or
> other protocol information is required for read access.

## Once access is granted

On the cluster, in the W1 environment:

```bash
cd "$BRP_PROJECT_ROOT/BetterRiseProject"
conda activate brp-w1
pip install -e '.[a2]'
synapse config          # paste a View-scope personal access token when prompted
chmod 600 ~/.synapseConfig
python -m src.reference.jobs.a2_synapse_inventory --no-write
```

**Authentication is a personal access token, not a password.** `synapse login
--rememberMe` was written here from older documentation and the flag does not
exist in synapseclient 4, which `pyproject.toml` pins. Make the token at
synapse.org under Account Settings, Personal Access Tokens, and give it **View
scope only** — a metadata inventory needs nothing else, and Download or Modify
scope on a shared cluster is a credential that can do more than the job it was
made for.

`synapse config` writes `~/.synapseConfig` in your home directory, outside the
repository. Keep it there. Do not export the token into a shell environment you
will leave running, and never place it anywhere under the working tree.

**The authenticated test is complete.** Anonymous refusal alone would not
establish an access gate, but the configured account received HTTP 403 for one
of the three root entities on 2026-09-09. Do not rerun the inventory until Sage
confirms access has changed; it cannot produce an honest complete inventory in
its current state.

The last command produces a metadata-only inventory. Review candidate products
for all four required properties before downloading anything:

1. stable cell IDs;
2. cell coordinates or polygons;
3. per-cell CDX2 intensity and the other panel markers; and
4. a documented segmentation/provenance field.

If all four exist, freeze the table identifier and write the A2 equivalence
pre-registration. If any is absent, record that explicitly and make a separate
segmentation decision; do not repurpose the public pixel arrays into a
post-hoc substitute.
