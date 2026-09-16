"""Bundle only the anonymous LaTeX source graph and referenced figure assets."""
from pathlib import Path
import json
import hashlib
from zipfile import ZipFile, ZIP_DEFLATED
from check_submission import ROOT, source_graph, GRAPHICS, main as check

def main():
    check()
    graph=source_graph(ROOT/'main.tex')
    files=set(graph+[ROOT/'refs.bib',ROOT/'neurips_2026.sty',ROOT/'main.bbl'])
    for path in graph:files.update(ROOT/name for name in GRAPHICS.findall(path.read_text()))
    files.update([ROOT/'diagnostics/low_count_diagnostic.csv',
                  ROOT/'diagnostics/candidate_reversals.json',
                  ROOT/'diagnostics/coverage_precision.json',
                  ROOT/'diagnostics/source_pool_properties.csv',
                  ROOT/'diagnostics/source_pool_summary.json',
                  ROOT/'evidence/review_record.json'])
    files.update(ROOT/f'diagnostics/external_control_{name}.csv'
                 for name in ['primary','shift','balance'])
    files.update(ROOT/f'diagnostics/{name}.csv' for name in
                 ['precision_followup','precision_seed_rates','precision_paired',
                  'retained_control_strata','retained_control_bulk'])
    files.update(ROOT/f'sections/{name}.tex' for name in
                 ['refdesign','learned_table','matrix_table'])
    target=ROOT/'source.zip'
    with ZipFile(target,'w',ZIP_DEFLATED) as z:
        for path in sorted(files):
            if path.suffix in {'.tex','.bib','.bbl','.json','.csv'}:
                content=path.read_text().replace(r'\_', '_').lower()
                for token in ['gate_memo_w2','harness_design_spec','betterriseproject','bodebosell','diegomardian','/users/']:
                    assert token not in content, f'Identifying string in review bundle: {path.name}'
            z.write(path,path.relative_to(ROOT))
    result={'archive':'source.zip','sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
            'files':[str(p.relative_to(ROOT)) for p in sorted(files)]}
    (ROOT/'source_bundle.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Wrote source.zip with',len(files),'files')
if __name__=='__main__':main()
