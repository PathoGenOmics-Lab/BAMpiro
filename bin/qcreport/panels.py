"""Optional analysis panels built from the parsed per-sample VCFs.

build_dynamics turns metadata + variants into per-connected-group allele-frequency
trajectories; build_epistasis scores pairs of those trajectories for co-variation;
build_snp_matrix flattens the same variants into the sparse sample x site matrix.
Each returns None when there is nothing to show, and the report hides the panel.

_pearson / _perm_p / _amp are the statistics build_epistasis is built on.
"""
from __future__ import annotations

import random
import re
import zlib

_DYN_NONSYN    = re.compile(r'missense|stop_gained|stop_lost|start_lost|frameshift|inframe|splice|initiator', re.I)


def build_dynamics(metadata, variants, sample_meta=None, emerge=0.25, fix=0.90, loss=0.10, min_points=2, min_move=0.15):
    """Per-connected-group AF trajectories over time, only for SNPs that move. None if nothing to show.
    sample_meta (parse_sample_meta output) attaches each group's metadata signature so the report can
    filter the series by any group-invariant field (e.g. patient / lineage)."""
    if not metadata:
        return None
    groups = {}
    for s, md in metadata.items():
        g = md.get('group')
        if g and s in variants:
            groups.setdefault(g, []).append(s)
    out_groups = []
    for g, samples in sorted(groups.items()):
        samples = sorted(samples, key=lambda s: (metadata[s]['tnum'] is None,
                                                 metadata[s]['tnum'] if metadata[s]['tnum'] is not None else 0,
                                                 str(metadata[s]['time'])))
        times = [metadata[s]['time'] for s in samples]
        if len(set(times)) < min_points:
            continue
        allpos = set()
        for s in samples:
            allpos.update(variants[s].keys())
        series = []
        for pos in sorted(allpos):   # deterministic order (a set iterates in a PYTHONHASHSEED-dependent order)
            traj, dps, meta = [], [], None
            for s in samples:
                v = variants[s].get(pos)
                traj.append(v['af'] if v else 0.0)
                dps.append(v.get('dp') if v else None)   # per-timepoint depth (None where not called)
                if v and meta is None:
                    meta = v
            if max(traj) - min(traj) < min_move:
                continue
            flags = []
            if traj[0] <= loss and max(traj) >= emerge:
                flags.append('emergence')
            if traj[-1] >= fix and traj[0] < fix:
                flags.append('fixation')
            if traj[0] >= emerge and traj[-1] <= loss:
                flags.append('loss')
            if meta and _DYN_NONSYN.search(meta.get('eff', '')):
                flags.append('nonsyn')
            if meta and meta.get('imp', '') == 'HIGH':
                flags.append('high_impact')
            series.append({'pos': pos, 'gene': (meta or {}).get('gene', ''), 'eff': (meta or {}).get('eff', ''),
                           'imp': (meta or {}).get('imp', ''), 'alt': (meta or {}).get('alt', ''),
                           'aa': (meta or {}).get('aa', ''), 'aa_h37rv': (meta or {}).get('aa_h37rv', ''),
                           'gene_h37rv': (meta or {}).get('gene_h37rv', ''),
                           'pos_h37rv': (meta or {}).get('pos_h37rv', ''),
                           'traj': [round(x, 4) for x in traj], 'dp': dps, 'flags': flags})
        if not series:
            continue
        series.sort(key=lambda x: (-(max(x['traj']) - min(x['traj'])), -len(x['flags'])))
        gmeta = {}
        if sample_meta and sample_meta.get('fields'):
            rows = sample_meta.get('rows', {})
            for f in sample_meta['fields']:
                vals = sorted({v for v in ((rows.get(s) or {}).get(f) for s in samples) if v})
                if vals:
                    gmeta[f] = vals
        out_groups.append({'group': g, 'samples': samples, 'times': times,
                           'tnums': [metadata[s]['tnum'] for s in samples], 'meta': gmeta,
                           'series': series, 'n_flagged': sum(1 for x in series if x['flags'])})
    if not out_groups:
        return None
    return {'groups': out_groups, 'thresholds': {'emerge': emerge, 'fix': fix, 'loss': loss}}


def _pearson(x, y):
    """Pearson correlation of two equal-length trajectories, or None if undefined (n<3 or a flat one)."""
    n = len(x)
    if n < 3 or len(y) != n:
        return None
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx <= 1e-9 or syy <= 1e-9:
        return None
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return sxy / ((sxx * syy) ** 0.5)


def _perm_p(series_traj, obs_mean, B=2000, seed=0):
    """Permutation p-value for a pair's mean correlation: shuffle one trajectory's values within EACH
    series independently (breaking the time linkage) and recompute the mean |r|; p = fraction of the B
    nulls reaching the observed |r| (+1 smoothing). Deterministic (fixed per-pair seed) -> reproducible.
    With few timepoints the resolution is inherently coarse (a 4-point series has only 4!=24 orders), so
    a lone short series lands near p~0.04 while a pattern recurring across series drives p far lower."""
    obs = abs(obs_mean)
    rng = random.Random(seed)
    ge = 0
    for _ in range(B):
        tot, cnt = 0.0, 0
        for ta, tb in series_traj:
            perm = list(ta)
            rng.shuffle(perm)
            r = _pearson(perm, tb)
            if r is not None:
                tot += r
                cnt += 1
        if cnt and abs(tot / cnt) >= obs:
            ge += 1
    return (ge + 1) / (B + 1)


def _amp(traj):
    """Amplitude (max-min) of a trajectory's non-null allele frequencies; 0 if flat/empty. Used to pick the
    most-variable variants when a series has too many to test all pairs of."""
    vals = [v for v in (traj or []) if v is not None]
    return (max(vals) - min(vals)) if len(vals) >= 2 else 0.0


def build_epistasis(dynamics, min_r=0.8, min_points=3, top=300, perm=2000, max_vars=250, max_pairs_perm=800):
    """Candidate epistatic / linked variant pairs: two SNPs whose allele-frequency trajectories co-vary
    WITHIN a connected series - concordant (rise/fall together) or discordant (one rises as the other
    falls). Score = mean Pearson r of the two trajectories across every series where both move and the
    series has >= min_points timepoints. Significance = a permutation p-value (per pair), then a
    Benjamini-Hochberg FDR q-value across all reported pairs. Recurrence across independent series is
    tracked so a pattern seen in several series (low q) is separable from a lucky single short series.
    None if nothing qualifies. Built from the dynamics payload (moving variants per series).

    Scaling: the pairwise step is O(variants^2) per series and the permutation test is the dominant cost,
    so two caps keep it bounded for large cohorts (many variants and/or many samples -> long trajectories):
    each series contributes at most max_vars variants (the highest-amplitude / most-variable ones), and only
    the max_pairs_perm strongest-|r| candidate pairs get the (expensive) permutation p-value. Both caps are
    reported in the payload so the UI can say what was covered; the full trajectories stay in the dynamics/TSV."""
    if not dynamics or not dynamics.get('groups'):
        return None
    pairs, n_series = {}, 0
    vars_capped, n_vars_max = False, 0
    for g in dynamics['groups']:
        times = g.get('times') or []
        if len(set(times)) < min_points:
            continue
        n_series += 1
        series = g.get('series') or []
        n_vars_max = max(n_vars_max, len(series))
        if len(series) > max_vars:   # too many variants to test every pair: keep the most-variable ones
            vars_capped = True
            series = sorted(series, key=lambda s: -_amp(s.get('traj')))[:max_vars]
        for i in range(len(series)):
            for j in range(i + 1, len(series)):
                a, b = series[i], series[j]
                r = _pearson(a['traj'], b['traj'])
                if r is None:
                    continue
                lo, hi = (a, b) if str(a['pos']) <= str(b['pos']) else (b, a)  # order-independent key
                key = (lo['pos'], hi['pos'])
                rec = pairs.get(key)
                if rec is None:
                    rec = pairs[key] = {'A': lo, 'B': hi, 'rs': [], 'straj': [], 'bestn': -1,
                                        'times': times, 'trajA': lo['traj'], 'trajB': hi['traj'], 'group': g['group']}
                rec['rs'].append(r)
                rec['straj'].append((lo['traj'], hi['traj']))   # every series' pair, for the permutation test
                if len(set(times)) > rec['bestn']:   # keep the richest series for the mini-chart
                    rec['bestn'] = len(set(times))
                    rec['times'], rec['trajA'], rec['trajB'], rec['group'] = times, lo['traj'], hi['traj'], g['group']
    # Collect the pairs that clear the |r| threshold, then permutation-test only the strongest ones: the
    # permutation p-value is by far the dominant cost, so with many candidate pairs (common when trajectories
    # are short and spurious high correlations abound) we cap it to the top max_pairs_perm by |mean r|.
    cands = []
    for key, rec in pairs.items():
        rs = rec['rs']
        mean_r = sum(rs) / len(rs)
        if abs(mean_r) < min_r:
            continue
        cands.append((abs(mean_r), key, rec, mean_r))
    n_candidates = len(cands)
    pairs_capped = n_candidates > max_pairs_perm
    if pairs_capped:
        cands.sort(key=lambda c: -c[0])
        cands = cands[:max_pairs_perm]
    out = []
    for _absr, key, rec, mean_r in cands:
        rs = rec['rs']
        A, B = rec['A'], rec['B']
        n_pos = sum(1 for r in rs if r > 0)
        # seed from the STABLE pair key (positions), not the iteration index, so the permutation p-value is
        # reproducible regardless of dict/set iteration order. crc32 is deterministic (unlike hash()).
        seed = zlib.crc32(('%s|%s' % (key[0], key[1])).encode('utf-8'))
        p = _perm_p(rec['straj'], mean_r, B=perm, seed=seed)
        out.append({'geneA': A.get('gene', ''), 'posA': A['pos'], 'aaA': A.get('aa', ''), 'aaA_h37rv': A.get('aa_h37rv', ''), 'effA': A.get('eff', ''),
                    'geneB': B.get('gene', ''), 'posB': B['pos'], 'aaB': B.get('aa', ''), 'aaB_h37rv': B.get('aa_h37rv', ''), 'effB': B.get('eff', ''),
                    'r': round(mean_r, 3), 'n': len(rs), 'rmin': round(min(rs), 3), 'rmax': round(max(rs), 3),
                    'direction': 'concordant' if mean_r > 0 else 'discordant',
                    'p': p, 'consistent': (n_pos == len(rs) or n_pos == 0),
                    'times': rec['times'], 'trajA': rec['trajA'], 'trajB': rec['trajB'], 'group': rec['group']})
    if not out:
        return None
    # Benjamini-Hochberg FDR across every reported pair (multiple-testing correction)
    m = len(out)
    order = sorted(range(m), key=lambda i: out[i]['p'])
    qmin = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        qmin = min(qmin, out[i]['p'] * m / rank)
        out[i]['q'] = round(min(qmin, 1.0), 4)
    for pr in out:
        pr['tier'] = 'strong' if pr['q'] <= 0.05 else ('moderate' if pr['p'] <= 0.05 else 'weak')
        pr['p'] = round(pr['p'], 4)
    out.sort(key=lambda p: (p['q'], -abs(p['r'])))
    out = out[:top]
    # all-vs-all correlation matrix over the most-connected variants (INCLUDES sub-threshold cells, so it
    # is the full co-dynamics landscape, not only the reported pairs). Capped to bound the embedded size.
    allmeta, deg, meanr = {}, {}, {}
    for (a, b), rec in pairs.items():
        meanr[(a, b)] = sum(rec['rs']) / len(rec['rs'])
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
        for v in (rec['A'], rec['B']):
            allmeta.setdefault(v['pos'], {'pos': v['pos'], 'gene': v.get('gene', ''), 'aa': v.get('aa', '')})
    MAXNODES = 50
    node_pos = sorted(deg, key=lambda x: (-deg[x], x))[:MAXNODES]          # keep the best-connected variants
    node_pos.sort(key=lambda x: (allmeta[x]['gene'] or 'zzz', x))         # display order: gene, then position
    nidx = {p: i for i, p in enumerate(node_pos)}
    cells = {}
    for (a, b), mr in meanr.items():
        if a in nidx and b in nidx:
            i, j = nidx[a], nidx[b]
            lo, hi = (i, j) if i <= j else (j, i)
            cells['%d,%d' % (lo, hi)] = round(mr, 3)
    matrix = {'nodes': [allmeta[p] for p in node_pos], 'cells': cells,
              'truncated': len(deg) > MAXNODES, 'total_nodes': len(deg)}
    return {'pairs': out, 'min_r': min_r, 'min_points': min_points, 'n_series': n_series, 'perm': perm,
            'matrix': matrix,
            'vars_capped': vars_capped, 'max_vars': max_vars, 'n_vars_max': n_vars_max,
            'pairs_capped': pairs_capped, 'n_candidates': n_candidates, 'max_pairs_perm': max_pairs_perm,
            'n_concordant': sum(1 for p in out if p['direction'] == 'concordant'),
            'n_discordant': sum(1 for p in out if p['direction'] == 'discordant'),
            'n_strong': sum(1 for p in out if p['tier'] == 'strong')}


def build_snp_matrix(variants, reference='', max_sites=50000):
    """Sparse SNP matrix for the report: samples + one row per SNP site with its ref/alt/annotation
    and only the cells (sample index -> [af, dp]) actually called. The payload is gzip-compressed in the
    HTML, so this cap is a high safety backstop (most-shared sites kept) rather than a size limit; the
    pipeline also writes the complete matrix TSV separately."""
    if not variants:
        return None
    samples = list(variants.keys())
    sidx = {s: i for i, s in enumerate(samples)}
    sites = {}
    for s, posmap in variants.items():
        for key, v in posmap.items():
            contig, _, pos = key.rpartition(':')
            try:
                ipos = int(pos)
            except ValueError:
                continue
            st = sites.setdefault(key, {'contig': contig, 'pos': ipos, 'ref': v.get('ref', ''),
                                        'alt': set(), 'gene': '', 'eff': '', 'aa': '', 'aa_h37rv': '',
                                        'gene_h37rv': '', 'pos_h37rv': '', 'cells': {}})
            if v.get('alt'):
                st['alt'].add(v['alt'])
            for k in ('gene', 'eff', 'aa', 'aa_h37rv', 'gene_h37rv', 'pos_h37rv'):
                if v.get(k) and not st[k]:
                    st[k] = v[k]
            st['cells'][sidx[s]] = [v['af'], v.get('dp')]
    ordered = sorted(sites.values(), key=lambda x: (x['contig'], x['pos']))
    truncated = len(ordered) > max_sites
    if truncated:
        ordered = sorted(sites.values(), key=lambda x: -len(x['cells']))[:max_sites]
        ordered.sort(key=lambda x: (x['contig'], x['pos']))
    rows = [{'contig': x['contig'], 'pos': x['pos'], 'ref': x['ref'], 'alt': ','.join(sorted(x['alt'])),
             'gene': x['gene'], 'eff': x['eff'], 'aa': x['aa'], 'aa_h37rv': x['aa_h37rv'],
             'gene_h37rv': x['gene_h37rv'], 'pos_h37rv': x['pos_h37rv'], 'n': len(x['cells']), 'cells': x['cells']}
            for x in ordered]
    return {'samples': samples, 'reference': reference, 'rows': rows,
            'total_sites': len(sites), 'truncated': truncated}
