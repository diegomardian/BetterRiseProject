# A2: request access, then inventory metadata only

This is the next A2 action after the public-product inventory. It seeks only
to establish whether the original Chen/HTAN release contains a MANDO-derived
cell table or segmentation. It does not request a scientific result and does
not authorize downloading raw images.

## Access request text

Submit through the Sage service desk linked from the [existing request for these
same entities](https://www.synapse.org/Synapse%3Asyn5637528/discussion/threadId%3D11828).
Replace the bracketed fields before sending:

> **Subject:** Read-only metadata access request: COLON MAP MxIF / HTA11
>
> I am [name, institution, role] reusing the published COLON MAP / Chen et al.
> HTAN data for an academic, non-commercial reproducibility analysis. Please
> grant read access to `syn23564801`, `syn23630431`, and `syn23520239`.
>
> My immediate use is metadata inventory only: I will enumerate entity names,
> file types, and annotations to determine whether MANDO-derived cell-level
> intensity tables or segmentation masks exist. I will not download raw images
> or redistribute data. If a suitable processed table exists, I will use it
> only for the described same-cohort protein validation, with patient-level
> aggregation and no attempt to infer a per-cell GUCA2A result from this panel.
>
> Please let me know if a data-use agreement, training, or additional protocol
> information is required.

## Once access is granted

On the cluster, in the W1 environment:

```bash
cd "$BRP_PROJECT_ROOT/BetterRiseProject"
conda activate brp-w1
pip install -e '.[a2]'
synapse config          # paste a personal access token when prompted
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

**Try this before submitting the access request above.** The inventory in
`docs/a2_mxif_inventory.md` records that *anonymous* read fails and that no
login was ever attempted from this machine — which does not establish that the
entities are access-gated. Many Synapse entities need only an authenticated
account and accepted terms. If the authenticated inventory returns rows, no
request is needed; if it returns 403, send the request that day.

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
