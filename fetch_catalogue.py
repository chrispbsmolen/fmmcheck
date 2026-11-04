"""Download the catalogue's front page and every scheme file it lists.

  python3 fetch_catalogue.py            index.html, data/entries.json, data/tensors/*.mpl.bz2, data/SNAPSHOT_DATE
  python3 fetch_catalogue.py --pages    also data/page_ranks.json (rank and construction stated on each format page)

Resumable (existing files are kept).  The files stay on your machine; this repository does not redistribute them."""
import re, html, json, os, sys, time, datetime, concurrent.futures, urllib.request
BASE = 'https://fmm.univ-lille.fr/'
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data'); os.makedirs(D + '/tensors', exist_ok=True)
def get(url, tries=4):
    for a in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r: return r.read()
        except Exception as ex:
            err = ex; time.sleep(2 * (a + 1))
    raise err
if not os.path.exists(D + '/index.html'):
    open(D + '/index.html', 'wb').write(get(BASE + 'index.html'))
    open(D + '/SNAPSHOT_DATE', 'w').write(datetime.date.today().isoformat() + '\n')
s = open(D + '/index.html', encoding='utf-8').read()
entries = {}
for row in re.findall(r'<tr[^>]*>(.*?)</tr>', s, re.S):
    m = re.search(r'(\d+)x(\d+)x(\d+)\.html', row)
    if not m: continue
    txt = html.unescape(re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', ' ', row))).strip()
    mm = re.search(r'⟨(\d+)×(\d+)×(\d+):(\d+)⟩\s+(\d+)\s+(\d+)\s+(.*?)\s+(\.\d+|\d+\.\d*)\s+(\d+\.\d+)\s+(.*)$', txt)
    fmt = '%sx%sx%s' % m.groups()
    if mm:
        entries[fmt] = dict(fmt=fmt, X=int(mm.group(1)), Y=int(mm.group(2)), Z=int(mm.group(3)), rank=int(mm.group(4)),
                            naive=int(mm.group(6)), source=re.sub(r'\s*\[\s*DOI\s*\]\s*$', '', mm.group(7)).strip())
json.dump(entries, open(D + '/entries.json', 'w'), indent=0)
print('formats in the front-page table:', len(entries), flush=True)
def fetch(fmt):
    out = D + '/tensors/%s_tensor.mpl.bz2' % fmt
    if os.path.exists(out): return fmt, 'cached'
    data = get(BASE + '%s_tensor.mpl.bz2' % fmt)
    if data[:3] != b'BZh': return fmt, 'not bzip2'
    open(out + '.part', 'wb').write(data); os.replace(out + '.part', out); return fmt, 'ok'
with concurrent.futures.ThreadPoolExecutor(4) as ex:
    res = list(ex.map(fetch, sorted(entries)))
print('scheme files:', sum(1 for _, r in res if r in ('ok', 'cached')), 'problems:', [x for x in res if x[1] not in ('ok', 'cached')])
if '--pages' in sys.argv:
    P = json.load(open(D + '/page_ranks.json')) if os.path.exists(D + '/page_ranks.json') else {}
    for i, fmt in enumerate(sorted(entries)):
        if fmt in P: continue
        try:
            t = html.unescape(re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', ' ', get(BASE + fmt + '.html').decode('utf-8', 'replace'))))
            m = re.search(r'Description of fast matrix multiplication algorithm: ⟨(\d+)×(\d+)×(\d+):(\d+)⟩', t)
            d = re.search(r'Algorithm definition (.*?) Algorithm description', t)
            P[fmt] = dict(page_rank=int(m.group(4)) if m else None, definition=d.group(1).strip()[:300] if d else None)
        except Exception as ex:
            P[fmt] = dict(page_rank=None, error=repr(ex)[:100])
        if (i + 1) % 100 == 0: json.dump(P, open(D + '/page_ranks.json', 'w'), indent=0); print('pages', i + 1, flush=True)
    json.dump(P, open(D + '/page_ranks.json', 'w'), indent=0)
