#!/usr/bin/env python3
"""Gray codes and de Bruijn sequences through three small alphabets, and
four laser-cut pieces that use them.

    python3 walks.py              # check everything, write cut-files/ and index.html
    python3 walks.py --kerf=0.18  # cut files for another kerf (default 0.13)
    python3 walks.py --check      # check only, write nothing

A Gray code lists every symbol of an alphabet once, each differing from the
last by one change. A de Bruijn sequence is one cyclic string in which every
word of n symbols turns up exactly once, so a window n wide tells you where
you are. The alphabets:

    layer memory    sites A B C of close-packed layers (../atomic-layer-memory)
    Logic Alphabet  16 truth tables of two variables, each a 2 x 2 square of
                    stems (Shea Zellweger's script, ../Logic-Research)
    grid script     28 letters of two grid steps (../gernreich-writing)

Cut files (P2S layers as in the staircase ruler):
    green  #008000  engrave (filled)
    blue   #0000ff  score
    orange #ff8000  cut, inner (first)
    black  #000000  cut, outline
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
GREEN, BLUE, ORANGE, BLACK = '#008000', '#0000ff', '#ff8000', '#000000'
BED_W, BED_H = 600, 305
KERF = 0.13


# ---------------------------------------------------------------- the walks

def de_bruijn(k, n):
    """The least de Bruijn sequence B(k, n), cyclic, as a list of 0..k-1
    (Fredricksen, Kessler and Maiorana: the Lyndon words in order)."""
    a, seq = [0] * (k * n + 1), []

    def db(t, p):
        if t > n:
            if n % p == 0:
                seq.extend(a[1:p + 1])
        else:
            a[t] = a[t - p]
            db(t + 1, p)
            for j in range(a[t - p] + 1, k):
                a[t] = j
                db(t + 1, t)
    db(1, 1)
    return seq


def gray(n):
    """The reflected binary Gray code on n bits."""
    return [i ^ (i >> 1) for i in range(1 << n)]


def windows(seq, n, cyclic=True):
    m = len(seq) if cyclic else len(seq) - n + 1
    s = seq + seq[:n - 1] if cyclic else seq
    return [tuple(s[i:i + n]) for i in range(m)]


def bitstr(v, n):
    return format(v, f'0{n}b')


# ------------------------------------------------- 1. layer memory: A B C

SITES = 'ABC'
UP = {'A': 'B', 'B': 'C', 'C': 'A'}          # one turn of the stacking


def spacer_stack(bits):
    """The atomic-layer-memory code: 0 = A, 1 = C, a B between every bit."""
    return 'B'.join('C' if b == '1' else 'A' for b in bits)


def moved(s, t):
    return sum(x != y for x, y in zip(s, t))


def counting(n, codes):
    """Stacks for counting 0..2^n-1 in the given code order, and how many
    layers move at each step."""
    st = [spacer_stack(bitstr(c, n)) for c in codes]
    return st, [moved(st[i], st[i + 1]) for i in range(len(st) - 1)]


def stack_de_bruijn(n):
    """A cyclic stack of A B C layers, no two neighbours alike (a close-packed
    layer cannot sit on the site below it), in which every possible window of
    n layers turns up once. Such windows number 3 * 2^(n-1): an Euler circuit
    on the graph of (n-1)-layer windows (Hierholzer's walk)."""
    nodes = [w for w in _words(n - 1)]
    out = {v: [v[1:] + (s,) for s in SITES if s != v[-1]] for v in nodes}
    start = nodes[0]
    stack, circuit, it = [start], [], {v: 0 for v in nodes}
    while stack:
        v = stack[-1]
        if it[v] < len(out[v]):
            stack.append(out[v][it[v]])
            it[v] += 1
        else:
            circuit.append(stack.pop())
    circuit.reverse()
    return [v[-1] for v in circuit[1:]]


def _words(m):
    res = [(s,) for s in SITES]
    for _ in range(m - 1):
        res = [w + (s,) for w in res for s in SITES if s != w[-1]]
    return res


def steps(stack):
    """Read a stack the other way: each layer is one turn up (A->B->C->A),
    1, or one turn down, 0, from the layer below."""
    return ''.join('1' if UP[a] == b else '0' for a, b in zip(stack, stack[1:]))


# ------------------------------------------------- 2. the Logic Alphabet

ROWS = ('TT', 'TF', 'FT', 'FF')
MEANING = {
    '1111': 'always true', '0000': 'never true',
    '1000': 'A and B', '0111': 'not both A and B',
    '1110': 'A or B', '0001': 'neither A nor B',
    '1011': 'A implies B', '0100': 'A and not B',
    '1101': 'B implies A', '0010': 'B and not A',
    '1100': 'A', '0011': 'not A', '1010': 'B', '0101': 'not B',
    '1001': 'A if and only if B', '0110': 'A or B but not both',
}


def logic_cycles():
    """Every Gray cycle through the 16 truth tables (a Hamiltonian cycle on
    the 4-cube), as directed cycles from 0000."""
    res = []

    def dfs(path, seen):
        v = path[-1]
        if len(path) == 16:
            if bin(v).count('1') == 1:
                res.append(path[:])
            return
        for b in (0, 1, 2, 3):
            w = v ^ (1 << b)
            if not seen >> w & 1:
                path.append(w)
                dfs(path, seen | 1 << w)
                path.pop()
    dfs([0], 1)
    return res


def antipodal(c):
    return all(c[(i + 8) % 16] == c[i] ^ 15 for i in range(16))


def logic_ring():
    """The first antipodal Gray cycle found, turned to start at 1111: each
    truth table has its negation directly across the ring."""
    c = next(c for c in logic_cycles() if antipodal(c))
    i = c.index(15)
    return [bitstr(v, 4) for v in c[i:] + c[:i]]


# ------------------------------------------------- 3. the grid script

DIRS = ['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE']          # 45 deg apart
STEP = {'E': (1, 0), 'W': (-1, 0), 'N': (0, -1), 'S': (0, 1),
        'NE': (1, -1), 'NW': (-1, -1), 'SE': (1, 1), 'SW': (-1, 1)}
# the 28 letters as gernreich-writing/gernreich.js traces them
LETTERS = {
    'KH': ('N', 'SE'), 'Ř': ('E', 'SW'), 'Z': ('SE', 'N'), 'Q': ('SW', 'E'),
    'J': ('S', 'NE'), 'X': ('W', 'SE'), 'K': ('NE', 'S'), 'V': ('SE', 'W'),
    'B': ('S', 'E'), 'W': ('N', 'E'), 'P': ('E', 'S'), 'Y': ('S', 'W'),
    'G': ('E', 'E'), 'F': ('S', 'S'), 'M': ('SE', 'SW'), 'C': ('SE', 'NE'),
    'U': ('NE', 'SE'), 'L': ('SW', 'SE'), 'D': ('SW', 'S'), 'H': ('E', 'SE'),
    'S': ('S', 'SW'), 'I': ('SE', 'E'), 'R': ('SE', 'S'), 'N': ('E', 'NE'),
    'O': ('S', 'SE'), 'A': ('NE', 'E'), 'T': ('SE', 'SE'), 'E': ('NE', 'NE'),
}


def opp(d):
    return DIRS[(DIRS.index(d) + 4) % 8]


def eighths(a, b):
    """How far b is turned from a, in eighths of a turn, 0..4."""
    t = (DIRS.index(b) - DIRS.index(a)) % 8
    return min(t, 8 - t)


def bend(a, b):
    """Signed turn at the letter's middle, in eighths, -3..3."""
    return (DIRS.index(b) - DIRS.index(a) + 4) % 8 - 4


def tracings(p):
    a, b = p
    return [(a, b), (opp(b), opp(a))]          # a letter reads either way


SHAPE = {t: k for k, p in LETTERS.items() for t in tracings(p)}


def grid_moves():
    """Letters one change apart: some tracing of each keeps one step and
    turns the other. The cost is that turn in eighths."""
    nb = {k: {} for k in LETTERS}
    for k, p in LETTERS.items():
        for a, b in tracings(p):
            for c in DIRS:
                for q, cost in (((c, b), eighths(a, c)), ((a, c), eighths(b, c))):
                    o = SHAPE.get(q)
                    if o and o != k:
                        nb[k][o] = min(nb[k].get(o, 9), cost)
    return nb


def grid_sides():
    """Turning one step by an eighth changes the bend by one, so letters of
    odd bend and even bend alternate along such moves."""
    odd = sorted(k for k, (a, b) in LETTERS.items() if bend(a, b) % 2)
    even = sorted(k for k in LETTERS if k not in odd)
    return even, odd


def grid_bound():
    """A cycle through all 28 crosses between the sides an even number of
    times, at most twice the smaller side; every other move stays on its
    side and costs at least two eighths."""
    even, odd = grid_sides()
    cross = 2 * min(len(even), len(odd))
    return cross + 2 * (len(LETTERS) - cross)


def grid_ring():
    """The first Gray cycle through the 28 letters, from G, that meets
    grid_bound(): the least turning possible."""
    nb, keys, bound = grid_moves(), list(LETTERS), grid_bound()
    found = []

    def dfs(path, seen, cost):
        if found or cost + (len(keys) - len(path)) > bound:
            return
        v = path[-1]
        if len(path) == len(keys):
            if path[0] in nb[v] and cost + nb[v][path[0]] == bound:
                found.append(path[:])
            return
        for w, c in sorted(nb[v].items(), key=lambda x: (x[1], x[0])):
            if w not in seen:
                path.append(w)
                seen.add(w)
                dfs(path, seen, cost + c)
                seen.discard(w)
                path.pop()
    dfs(['G'], {'G'}, 0)
    return found[0]


def grid_change(k, o):
    """How letter k becomes o: the tracing of k, the tracing of o, the slot
    turned (0 or 1) and the turn in eighths."""
    best = None
    for a, b in tracings(LETTERS[k]):
        for q in tracings(LETTERS[o]):
            for slot in (0, 1):
                if (a, b)[1 - slot] == q[1 - slot]:
                    t = eighths((a, b)[slot], q[slot])
                    if best is None or t < best[3]:
                        best = ((a, b), q, slot, t)
    return best


# ---------------------------------------------------------------- checks

def check():
    out = []

    def ok(cond, what):
        out.append(('PASS ' if cond else 'FAIL ') + what)
        return cond

    for k, n in ((2, 5), (2, 9)):
        s = de_bruijn(k, n)
        w = windows(s, n)
        ok(len(s) == k ** n and len(set(w)) == k ** n,
           f'B({k},{n}) has length {k ** n} and every window of {n} once')
    track = ruler_track(300)
    ok(len(set(windows(track, 9, cyclic=False))) == 301,
       'the ruler track: 301 windows of 9 cells, all different')
    for n in (5,):
        g = gray(n)
        ok(all(bin(g[i] ^ g[(i + 1) % len(g)]).count('1') == 1 for i in range(len(g))),
           f'the {n}-bit Gray code changes one bit per step, round the circle')
    st, mv = counting(4, range(16))
    gst, gmv = counting(4, gray(4))
    ok(sum(mv) == 26 and sum(gmv) == 15 and max(gmv) == 1,
       f'counting 0..15 in the layer stack moves {sum(mv)} layers in binary, {sum(gmv)} in Gray')
    for n in (3, 4, 5):
        s = stack_de_bruijn(n)
        w = windows(s, n)
        good = all(a != b for a, b in zip(s, s[1:] + s[:1]))
        ok(good and len(s) == 3 * 2 ** (n - 1) and len(set(w)) == len(s) == len(_words(n)),
           f'layer stack, windows of {n}: {len(s)} layers, no two neighbours alike, every window once')
    s = stack_de_bruijn(4)
    ok(steps(s) and all(x != y for x, y in zip(s, s[1:])), 'steps() reads every layer after the first')
    cyc = logic_cycles()
    ok(len(cyc) == 2688, f'{len(cyc) // 2} Gray cycles through the 16 truth tables')
    anti = [c for c in cyc if antipodal(c)]
    ok(len(anti) == 96, f'{len(anti) // 2} of them antipodal')
    ring = logic_ring()
    ok(sorted(ring) == sorted(MEANING) and
       all(moved(ring[i], ring[(i + 1) % 16]) == 1 for i in range(16)) and
       all(int(ring[i], 2) ^ int(ring[i + 8], 2) == 15 for i in range(8)),
       'the Logic ring: all 16, one stem per step, each across from its negation')
    ok(len(SHAPE) == 56 and len(set(SHAPE.values())) == 28, '28 letters, 56 tracings')
    even, odd = grid_sides()
    ok((len(even), len(odd)) == (12, 16), f'grid letters: {len(even)} of even bend, {len(odd)} of odd')
    nb = grid_moves()
    one = {k: {o for o, c in v.items() if c == 1} for k, v in nb.items()}
    ok(all((k in even) != (o in even) for k in one for o in one[k]),
       'every one-eighth move crosses between the sides, so no cycle uses only those')
    gr = grid_ring()
    cost = [nb[gr[i]][gr[(i + 1) % 28]] for i in range(28)]
    ok(sorted(gr) == sorted(LETTERS) and sum(cost) == grid_bound() == 32,
       f'the grid ring: all 28, {cost.count(1)} one-eighth turns and {cost.count(2)} quarter turns, '
       f'{sum(cost)} eighths, the least possible')
    for line in out:
        print(line)
    return all(l.startswith('PASS') for l in out)


# ------------------------------------------------------------ cut files

def f(v):
    s = f'{v:.3f}'.rstrip('0').rstrip('.')
    return '0' if s == '-0' else s


DIGITS = {  # single-stroke digits on a 1 x 2 box, y down (from staircase_ruler.py)
    '0': [[(.2, 0), (.8, 0), (1, .3), (1, 1.7), (.8, 2), (.2, 2), (0, 1.7), (0, .3), (.2, 0)]],
    '1': [[(.25, .4), (.6, 0), (.6, 2)]],
    '2': [[(0, .4), (.25, 0), (.75, 0), (1, .4), (1, .8), (0, 2), (1, 2)]],
    '3': [[(0, 0), (1, 0), (.4, .8), (.75, .8), (1, 1.1), (1, 1.7), (.75, 2), (.25, 2), (0, 1.7)]],
    '4': [[(.75, 2), (.75, 0), (0, 1.4), (1, 1.4)]],
    '5': [[(1, 0), (.1, 0), (0, .9), (.7, .8), (1, 1.1), (1, 1.7), (.75, 2), (.25, 2), (0, 1.7)]],
    '6': [[(.9, 0), (.4, 0), (0, .6), (0, 1.7), (.25, 2), (.75, 2), (1, 1.7), (1, 1.2), (.75, .9),
           (.25, .9), (0, 1.2)]],
    '7': [[(0, 0), (1, 0), (.35, 2)]],
    '8': [[(.5, .9), (.1, .7), (.1, .25), (.35, 0), (.65, 0), (.9, .25), (.9, .7), (.5, .9),
           (0, 1.2), (0, 1.7), (.25, 2), (.75, 2), (1, 1.7), (1, 1.2), (.5, .9)]],
    '9': [[(.1, 2), (.6, 2), (1, 1.4), (1, .3), (.75, 0), (.25, 0), (0, .3), (0, .8), (.25, 1.1),
           (.75, 1.1), (1, .8)]],
}


def number(text, cx, cy, h):
    w, gap = h / 2, h * 0.15
    x0 = cx - (len(text) * w + (len(text) - 1) * gap) / 2
    return [[(x0 + i * (w + gap) + px * w, cy - h / 2 + py * w) for px, py in s]
            for i, ch in enumerate(text) for s in DIGITS[ch]]


def pl(pts, close=False):
    return 'M' + ' L'.join(f'{f(x)},{f(y)}' for x, y in pts) + (' Z' if close else '')


def circle(cx, cy, r):
    return (f'M{f(cx - r)},{f(cy)} A{f(r)},{f(r)} 0 1 0 {f(cx + r)},{f(cy)} '
            f'A{f(r)},{f(r)} 0 1 0 {f(cx - r)},{f(cy)} Z')


def rect(x, y, w, h):
    return f'M{f(x)},{f(y)} H{f(x + w)} V{f(y + h)} H{f(x)} Z'


def polar(cx, cy, r, turn):
    """Point at radius r, `turn` of a full turn clockwise from the top."""
    a = 2 * math.pi * turn
    return cx + r * math.sin(a), cy - r * math.cos(a)


def wedge(cx, cy, r0, r1, t0, t1):
    """Annular sector from turn t0 to t1 clockwise, radii r0 < r1."""
    big = 1 if t1 - t0 > 0.5 else 0
    a, b = polar(cx, cy, r1, t0), polar(cx, cy, r1, t1)
    c, d = polar(cx, cy, r0, t1), polar(cx, cy, r0, t0)
    return (f'M{f(a[0])},{f(a[1])} A{f(r1)},{f(r1)} 0 {big} 1 {f(b[0])},{f(b[1])} '
            f'L{f(c[0])},{f(c[1])} A{f(r0)},{f(r0)} 0 {big} 0 {f(d[0])},{f(d[1])} Z')


def runs(cells):
    """(start, length) of each run of 1s in a list."""
    res, i = [], 0
    while i < len(cells):
        if cells[i]:
            j = i
            while j < len(cells) and cells[j]:
                j += 1
            res.append((i, j - i))
            i = j
        else:
            i += 1
    return res


class Piece:
    """Paths by layer; written as a cut file or as a page preview."""

    def __init__(self, name, box):
        self.name, self.box = name, box            # box: x0, y0, x1, y1
        self.engrave, self.score, self.inner, self.outline = [], [], [], []

    def size(self):
        x0, y0, x1, y1 = self.box
        return x1 - x0, y1 - y0

    def cut_svg(self):
        x0, y0, x1, y1 = self.box
        p = 1
        W, H = x1 - x0 + 2 * p, y1 - y0 + 2 * p
        g = lambda col, ds: (f'  <g fill="none" stroke="{col}" stroke-width="0.1">\n'
                             + ''.join(f'    <path d="{d}"/>\n' for d in ds) + '  </g>\n')
        return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{f(W)}mm" height="{f(H)}mm" '
                f'viewBox="{f(x0 - p)} {f(y0 - p)} {f(W)} {f(H)}">\n'
                f'  <g fill="{GREEN}" stroke="none">\n'
                + ''.join(f'    <path d="{d}"/>\n' for d in self.engrave) + '  </g>\n'
                + g(BLUE, self.score) + g(ORANGE, self.inner) + g(BLACK, self.outline)
                + '</svg>\n')

    def preview(self, label, extra=''):
        x0, y0, x1, y1 = self.box
        p = 1
        W, H = x1 - x0 + 2 * p, y1 - y0 + 2 * p
        body = ' '.join(self.outline + self.inner)
        return (f'<svg class="piece" viewBox="{f(x0 - p)} {f(y0 - p)} {f(W)} {f(H)}" role="img" '
                f'aria-label="{label}">'
                f'<path class="stock" fill-rule="evenodd" d="{body}"/>'
                + ''.join(f'<path class="eng" d="{d}"/>' for d in self.engrave)
                + f'<path class="sco" d="{" ".join(self.score)}"/>' + extra + '</svg>')


# 1. the window ruler -------------------------------------------------------

RULER_W, TRACK_H, TAIL, HOLE_D, CORNER_R, BLEED = 30, 8, 20, 6, 4, 0.3
WINDOW = 9


def ruler_track(length):
    """Cells of the ruler's track, 1mm each: the start of B(2,9) unrolled,
    long enough that a 9mm window at every whole mm 0..length sits on it."""
    s = de_bruijn(2, WINDOW)
    return (s + s[:WINDOW - 1])[:length + WINDOW]


def window_ruler(length, kerf):
    k = kerf / 2
    W = RULER_W
    pc = Piece(f'window-ruler-{length}', (-k, -k, length + TAIL + k, W + k))
    track = ruler_track(length)
    ty = W - TRACK_H
    for x, n in runs(track):
        pc.engrave.append(rect(x, ty, n, TRACK_H + BLEED))
    for mm in range(0, len(track) + 1):
        pc.score.append(pl([(mm, ty), (mm, ty - (3 if mm % 10 == 0 else 1.5))]))
    for mm in range(1, length + 1):
        n = mm % 10
        pc.score.append(pl([(mm, 0), (mm, 7 if n == 0 else 5 if n == 5 else 3)]))
    for mm in range(10, length + 1, 10):
        pc.score += [pl(s) for s in number(str(mm), mm, 11.5, 4)]
    x0, y0, x1, y1, r = -k, -k, length + TAIL + k, W + k, CORNER_R + k
    pc.outline.append(f'M{f(x0)},{f(y0)} H{f(x1 - r)} A{f(r)},{f(r)} 0 0 1 {f(x1)},{f(y0 + r)} '
                      f'V{f(y1 - r)} A{f(r)},{f(r)} 0 0 1 {f(x1 - r)},{f(y1)} H{f(x0)} Z')
    pc.inner.append(circle(length + TAIL / 2, W / 2, HOLE_D / 2 - k))
    return pc


def window_reader(kerf):
    """A card to lay on the track: a window 9 cells wide, ticks at every
    cell edge, and a line from the window's left edge to read the scale by."""
    k = kerf / 2
    w, h = 36, 30
    pc = Piece('window-reader', (-k, -k, w + k, h + k))
    wx, wy = (w - WINDOW) / 2, h - 4 - TRACK_H
    pc.inner.append(rect(wx + k, wy + k, WINDOW - 2 * k, TRACK_H - 2 * k))
    for i in range(WINDOW + 1):
        pc.score.append(pl([(wx + i, wy - 0.6), (wx + i, wy - (2.5 if i in (0, WINDOW) else 1.5))]))
    pc.score.append(pl([(wx, wy - 3), (wx, 0)]))
    pc.score.append(pl([(wx, wy + TRACK_H + 0.6), (wx, h)]))
    pc.outline.append(rect(-k, -k, w + 2 * k, h + 2 * k))
    return pc


# 2. two encoder discs --------------------------------------------------------

DISC_R, DISC_HOLE = 45, 6


def gray_disc(kerf, cx=0, cy=0):
    k = kerf / 2
    pc = Piece('gray-disc', (cx - DISC_R - k, cy - DISC_R - k, cx + DISC_R + k, cy + DISC_R + k))
    g, n = gray(5), 32
    for ring in range(5):
        r0, r1 = 12 + 6 * ring, 18 + 6 * ring
        cells = [(g[i] >> (4 - ring)) & 1 for i in range(n)]
        for s, m in runs(cells):
            pc.engrave.append(wedge(cx, cy, r0, r1, s / n, (s + m) / n))
        pc.score.append(circle(cx, cy, r1) if ring == 4 else circle(cx, cy, r0))
    for i in range(n):
        a, b = polar(cx, cy, 42.6, i / n), polar(cx, cy, 44 if i % 8 else 44.4, i / n)
        pc.score.append(pl([a, b]))
    pc.inner.append(circle(cx, cy, DISC_HOLE / 2 - k))
    pc.outline.append(circle(cx, cy, DISC_R + k))
    return pc


def de_bruijn_disc(kerf, cx=0, cy=0):
    k = kerf / 2
    pc = Piece('de-bruijn-disc', (cx - DISC_R - k, cy - DISC_R - k, cx + DISC_R + k, cy + DISC_R + k))
    s, n = de_bruijn(2, 5), 32
    for st, m in runs(s):
        pc.engrave.append(wedge(cx, cy, 30, 42, st / n, (st + m) / n))
    pc.score += [circle(cx, cy, 30), circle(cx, cy, 42)]
    for i in range(n):
        a, b = polar(cx, cy, 42.6, i / n), polar(cx, cy, 44 if i % 8 else 44.4, i / n)
        pc.score.append(pl([a, b]))
        pc.score.append(pl([polar(cx, cy, 28, i / n), polar(cx, cy, 29.4, i / n)]))
    pc.inner.append(circle(cx, cy, DISC_HOLE / 2 - k))
    pc.outline.append(circle(cx, cy, DISC_R + k))
    return pc


# 3. the Logic ring and 4. the grid ring --------------------------------------

LR_R, LR_r, LR_C, CELL, GAP = 62, 40, 51, 5, 1


def logic_ring_piece(kerf):
    k = kerf / 2
    pc = Piece('logic-ring', (-LR_R - k, -LR_R - k, LR_R + k, LR_R + k))
    for i, t in enumerate(logic_ring()):
        cx, cy = polar(0, 0, LR_C, i / 16)
        for j, b in enumerate(t):
            x = cx - CELL - GAP / 2 + (j % 2) * (CELL + GAP)
            y = cy - CELL - GAP / 2 + (j // 2) * (CELL + GAP)
            pc.score.append(rect(x, y, CELL, CELL))
            if b == '1':
                pc.engrave.append(rect(x + 0.4, y + 0.4, CELL - 0.8, CELL - 0.8))
    # this way up: the glyphs only mean what they mean upright
    pc.score.append(pl([(-2, -LR_r - 3.5), (0, -LR_r - 6), (2, -LR_r - 3.5)]))
    pc.inner.append(circle(0, 0, LR_r - k))
    pc.outline.append(circle(0, 0, LR_R + k))
    return pc


GR_R, GR_r, GR_C, GSTEP = 68, 44, 56, 2.6


def letter_points(trace, x=0, y=0, s=GSTEP):
    pts = [(x, y)]
    for d in trace:
        dx, dy = STEP[d]
        pts.append((pts[-1][0] + s * dx, pts[-1][1] + s * dy))
    return pts


def centred(trace, cx, cy, s=GSTEP):
    pts = letter_points(trace, 0, 0, s)
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    ox, oy = cx - (min(xs) + max(xs)) / 2, cy - (min(ys) + max(ys)) / 2
    return [(px + ox, py + oy) for px, py in pts]


def grid_ring_piece(kerf):
    k = kerf / 2
    pc = Piece('grid-ring', (-GR_R - k, -GR_R - k, GR_R + k, GR_R + k))
    for i, name in enumerate(grid_ring()):
        cx, cy = polar(0, 0, GR_C, i / 28)
        pts = centred(LETTERS[name], cx, cy)
        pc.score.append(pl(pts))
        pc.engrave.append(circle(pts[0][0], pts[0][1], 0.6))
    pc.score.append(pl([(-2, -GR_r - 3.5), (0, -GR_r - 6), (2, -GR_r - 3.5)]))
    pc.inner.append(circle(0, 0, GR_r - k))
    pc.outline.append(circle(0, 0, GR_R + k))
    return pc


def pieces(kerf):
    return [window_ruler(300, kerf), window_reader(kerf), gray_disc(kerf), de_bruijn_disc(kerf),
            logic_ring_piece(kerf), grid_ring_piece(kerf)]


# ---------------------------------------------------------------- the page

HI, TURNED = ' class="hi"', ' class="turned"'
COL = {'A': 'var(--a)', 'B': 'var(--b)', 'C': 'var(--c)'}


def stack_svg(st, hi=(), w=14, h=9, label=''):
    """One stack as a column, bottom layer at the bottom."""
    n = len(st)
    cells = ''.join(
        f'<rect x="1" y="{1 + (n - 1 - i) * h}" width="{w}" height="{h - 1}" rx="1.5" '
        f'fill="{COL[c]}"{HI if i in hi else ""}/>'
        f'<text x="{1 + w / 2}" y="{1 + (n - 1 - i) * h + h / 2 + 0.5}">{c}</text>'
        for i, c in enumerate(st))
    return (f'<svg class="stack" viewBox="0 0 {w + 2} {n * h + 1}" width="{w + 2}" '
            f'height="{n * h + 1}" role="img" aria-label="{label}">{cells}</svg>')


def counting_html(n=4):
    rows = []
    for name, codes in (('Binary', list(range(16))), ('Gray', gray(4))):
        st, mv = counting(n, codes)
        cols = []
        for i, s in enumerate(st):
            hi = () if i == 0 else tuple(j for j in range(len(s)) if s[j] != st[i - 1][j])
            cols.append(f'<div class="col">{stack_svg(s, hi, label=f"{bitstr(codes[i], n)}")}'
                        f'<code>{bitstr(codes[i], n)}</code>'
                        f'<span class="mv">{"" if i == 0 else mv[i - 1]}</span></div>')
        rows.append(f'<h3>{name}: {sum(mv)} layers moved</h3><div class="stacks">{"".join(cols)}</div>')
    return '\n'.join(rows)


def stack_db_html(n=3):
    s = stack_de_bruijn(n)
    cells = ''.join(f'<span class="L" style="background:{COL[c]}">{c}</span>' for c in s)
    return (f'<div class="row-stack" aria-label="a cyclic stack of {len(s)} layers">{cells}'
            f'<span class="L wrap">{s[0]}</span><span class="L wrap">{s[1]}</span></div>'), s


def logic_table():
    ring = logic_ring()
    tr = []
    for i, t in enumerate(ring):
        nxt = ring[(i + 1) % 16]
        row = next(ROWS[j] for j in range(4) if t[j] != nxt[j])
        sq = ''.join(f'<i class="{"on" if b == "1" else ""}"></i>' for b in t)
        tr.append(f'<tr><td class="num">{i}</td><td><span class="sq">{sq}</span></td>'
                  f'<td><code>{t}</code></td><td>{MEANING[t]}</td>'
                  f'<td>{"on" if nxt[ROWS.index(row)] == "1" else "off"} {row}</td>'
                  f'<td>{MEANING[ring[(i + 8) % 16]]}</td></tr>')
    return '\n'.join(tr)


def letter_svg(trace, hi_slot=None, s=10):
    pts = letter_points(trace, 0, 0, s)
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    x0, y0 = min(xs) - 4, min(ys) - 4
    w, h = max(xs) - min(xs) + 8, max(ys) - min(ys) + 8
    seg = ''.join(f'<line x1="{f(pts[i][0])}" y1="{f(pts[i][1])}" x2="{f(pts[i + 1][0])}" '
                  f'y2="{f(pts[i + 1][1])}"{TURNED if i == hi_slot else ""}/>'
                  for i in range(2))
    return (f'<svg class="letter" viewBox="{f(x0)} {f(y0)} {f(w)} {f(h)}" width="{f(w * 2.2)}" '
            f'height="{f(h * 2.2)}" aria-hidden="true">{seg}'
            f'<circle cx="{f(pts[0][0])}" cy="{f(pts[0][1])}" r="2"/></svg>')


def grid_list():
    ring = grid_ring()
    out = []
    for i, name in enumerate(ring):
        prev = ring[i - 1]
        a, q, slot, t = grid_change(prev, name)
        out.append(f'<div class="glet">{letter_svg(q, slot)}<b>{name}</b>'
                   f'<span>{"⅛" if t == 1 else "¼"}</span></div>')
    return ''.join(out)


def page(kerf):
    ps = {p.name: p for p in pieces(kerf)}
    track = ''.join(map(str, ruler_track(300)))
    disc = ''.join(map(str, de_bruijn(2, 5)))
    sdb, s3 = stack_db_html(3)
    even, odd = grid_sides()
    sub = {
        '{{RULER}}': ps['window-ruler-300'].preview('the window ruler'),
        '{{READER}}': ps['window-reader'].preview('the reader card'),
        '{{GRAY_DISC}}': ps['gray-disc'].preview('the Gray code disc'),
        '{{DB_DISC}}': ps['de-bruijn-disc'].preview('the de Bruijn disc'),
        '{{LOGIC_RING}}': ps['logic-ring'].preview('the Logic Alphabet ring'),
        '{{GRID_RING}}': ps['grid-ring'].preview('the grid-script ring'),
        '{{TRACK}}': track, '{{DISC}}': disc,
        '{{COUNTING}}': counting_html(),
        '{{STACK_DB}}': sdb,
        '{{STACK_DB_STEPS}}': steps(s3 + s3[:1]),
        '{{LOGIC_TABLE}}': logic_table(),
        '{{GRID_LIST}}': grid_list(),
        '{{EVEN}}': ' '.join(even), '{{ODD}}': ' '.join(odd),
        '{{LOGIC_JSON}}': json.dumps(logic_ring()),
        '{{STACK_JSON}}': json.dumps(s3),
        '{{MEANING_JSON}}': json.dumps(MEANING),
    }
    html = open(os.path.join(HERE, 'page.html')).read()
    for a, b in sub.items():
        assert a in html, a
        html = html.replace(a, b)
    assert '{{' not in html
    return html


def main(argv):
    kerf, write = KERF, True
    for a in argv:
        if a.startswith('--kerf='):
            kerf = float(a.split('=')[1])
        elif a == '--check':
            write = False
        else:
            sys.exit(f'unknown {a}; takes --kerf=MM, --check')
    if not check():
        sys.exit('checks failed')
    for p in pieces(kerf):
        w, h = p.size()
        fits = round(w, 3) <= BED_W and round(h, 3) <= BED_H
        print(f'{p.name:20s} {w:7.2f} x {h:6.2f} mm  {"fits" if fits else "DOES NOT FIT"} the P2S bed')
        if not fits:
            sys.exit(1)
        if write:
            os.makedirs(os.path.join(HERE, 'cut-files'), exist_ok=True)
            with open(os.path.join(HERE, 'cut-files', p.name + '.svg'), 'w') as fh:
                fh.write(p.cut_svg())
    if write:
        with open(os.path.join(HERE, 'index.html'), 'w') as fh:
            fh.write(page(KERF))
        print(f'wrote cut-files/ (kerf {kerf}) and index.html')


if __name__ == '__main__':
    main(sys.argv[1:])
