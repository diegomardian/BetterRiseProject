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
                  ROOT/'evidence/case_provenance.json'])
    target=ROOT/'source.zip'
    with ZipFile(target,'w',ZIP_DEFLATED) as z:
        for path in sorted(files):z.write(path,path.relative_to(ROOT))
    result={'archive':'source.zip','sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
            'files':[str(p.relative_to(ROOT)) for p in sorted(files)]}
    (ROOT/'source_bundle.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Wrote source.zip with',len(files),'files')
if __name__=='__main__':main()
