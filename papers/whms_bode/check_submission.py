"""Check the exact built PDF, its source graph, anonymity and official style."""
from pathlib import Path
import hashlib
import json
import re
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parent
INPUT=re.compile(r'\\(?:input|include)\{([^}]+)\}')
GRAPHICS=re.compile(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}')

def source_graph(path,seen=None):
    seen=set() if seen is None else seen
    path=path.resolve()
    if path in seen:return []
    seen.add(path)
    assert path.is_file(),f'Missing source: {path}'
    out=[path]
    for v in INPUT.findall(path.read_text()):
        child=ROOT/v
        out+=source_graph(child if child.suffix else child.with_suffix('.tex'),seen)
    return out

def main():
    graph=source_graph(ROOT/'main.tex')
    for path in graph:
        raw=path.read_text()
        assert not re.search(r'[\x00-\x08\x0b\x0c\x0e-\x1f]',raw),f'Control character in {path.name}'
        assert not re.search(r'(?<![A-Za-z])extbf',raw),f'Malformed bold command in {path.name}'
        for name in GRAPHICS.findall(raw):assert (ROOT/name).is_file(),name
    dependencies=graph+[ROOT/'refs.bib',ROOT/'neurips_2026.sty']
    for path in graph:
        dependencies.extend(ROOT/name for name in GRAPHICS.findall(path.read_text()))
    pdf_time=(ROOT/'main.pdf').stat().st_mtime
    assert all(path.stat().st_mtime <= pdf_time for path in dependencies),'PDF is stale; rebuild after source/figure edits'
    log=(ROOT/'main.log').read_text()
    bad=[line for line in log.splitlines() if line.startswith('!') or 'undefined' in line or 'Overfull' in line or 'Rerun to get' in line]
    assert not bad,'LaTeX defects:\n'+'\n'.join(bad)
    # The bibliography deliberately starts on a fresh page. Count every page
    # before its heading, so responsible-use spillover cannot pass unnoticed.
    reader=PdfReader(ROOT/'main.pdf')
    texts=[page.extract_text() or '' for page in reader.pages]
    ref=next(i for i,t in enumerate(texts) if re.match(r'^References(?:\d+)?\s',t))
    appendix=next(i for i,t in enumerate(texts) if re.match(r'^A\s+Simulation settings',t))
    assert appendix>ref,'Appendix must follow references'
    assert ref<=9,f'Main text is {ref} pages (limit 9)'
    assert 'Responsible use, limitations and impact' in '\n'.join(texts[:ref])
    assert 'Data and impact.' in '\n'.join(texts[:ref]),'Responsible-use statement is incomplete'
    combined='\n'.join(texts)+'\n'+'\n'.join(x.read_text() for x in graph)+'\n'+(ROOT/'refs.bib').read_text()
    for pat in [r'bodebosell',r'diegomardian',r'BetterRiseProject',r'/Users/',r'/home/[a-z]',r'extbfNone']:
        assert not re.search(pat,combined,re.I),f'Anonymity/typesetting issue: {pat}'
    meta=reader.metadata or {}
    assert not str(meta.get('/Author','')).strip(),'Author metadata present'
    assert not str(meta.get('/Keywords','')).strip(),'Keywords metadata present'
    source_meta=json.loads((ROOT/'style_provenance.json').read_text())
    digest=hashlib.sha256((ROOT/'neurips_2026.sty').read_bytes()).hexdigest()
    assert digest==source_meta['sha256'],'Official style has been modified'
    summary={'main_text_pages':ref,'reference_start_page':ref+1,
             'reference_pages':appendix-ref,'appendix_start_page':appendix+1,
             'appendix_pages':len(texts)-appendix,'total_pages':len(texts),
             'latex_errors':0,'undefined_references':0,'overfull_boxes':0,'anonymity':'passed',
             'official_style_sha256':digest,'pdf_sha256':hashlib.sha256((ROOT/'main.pdf').read_bytes()).hexdigest(),
             'source_files':[str(f.relative_to(ROOT)) for f in graph]}
    (ROOT/'build_verification.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(f"PASS: {ref} main-text pages; references start on {ref+1}; {len(texts)} total pages")
    print('PASS: no LaTeX errors, unresolved references, overfull boxes or detected identifying strings')
    print('PASS: official NeurIPS 2026 style checksum unchanged')
if __name__=='__main__':main()
