//! Gray codes and de Bruijn sequences through three small alphabets, and
//! four laser-cut pieces that use them.
//!
//!     target/release/walks              # check everything, write cut-files/ and index.html
//!     target/release/walks --kerf=0.18  # cut files for another kerf (default 0.13)
//!     target/release/walks --check      # check only, write nothing
//!     target/release/walks --resample   # rewrite grid-rings.json (seconds, not minutes)
//!     target/release/walks --graph      # the search graph, for tools/count_grid_rings.rs
//!
//! A Gray code lists every symbol of an alphabet once, each differing from the
//! last by one change. A de Bruijn sequence is one cyclic string in which every
//! word of n symbols turns up exactly once, so a window n wide tells you where
//! you are. The alphabets:
//!
//!     layer memory    sites A B C of close-packed layers (../atomic-layer-memory)
//!     Logic Alphabet  16 truth tables of two variables, each a 2 x 2 square of
//!                     stems (Shea Zellweger's script, ../Logic-Research)
//!     grid script     28 letters of two grid steps (../gernreich-writing)
//!
//! Cut files (P2S layers as in the staircase ruler):
//!     green  #008000  engrave (filled)
//!     blue   #0000ff  score
//!     orange #ff8000  cut, inner (first)
//!     black  #000000  cut, outline
//!
//! The Rust port of walks.py (2026-11-09): every file it writes, and
//! grid-rings.json from --resample, come out byte for byte as the Python's.
//! --resample draws from the same Mersenne Twister Python's random.Random(1)
//! does, so its shuffled search tries the neighbours in the same order.

use std::collections::{BTreeSet, HashMap, HashSet};
use std::f64::consts::PI;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

const GREEN: &str = "#008000";
const BLUE: &str = "#0000ff";
const ORANGE: &str = "#ff8000";
const BLACK: &str = "#000000";
const BED_W: f64 = 600.0;
const BED_H: f64 = 305.0;
/// how many least-turning grid rings exist in all
const GRID_COUNT: &str = "In all there are 2,161,152 least-turning rings, counted by tools/count_grid_rings.rs, a Rust search that took about 13 minutes.";
const KERF: f64 = 0.13;

fn here() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

// ---------------------------------------------------------------- the walks

/// The least de Bruijn sequence B(k, n), cyclic, as a list of 0..k-1
/// (Fredricksen, Kessler and Maiorana: the Lyndon words in order).
fn de_bruijn(k: usize, n: usize) -> Vec<usize> {
    fn db(t: usize, p: usize, k: usize, n: usize, a: &mut Vec<usize>, seq: &mut Vec<usize>) {
        if t > n {
            if n % p == 0 {
                seq.extend_from_slice(&a[1..p + 1]);
            }
        } else {
            a[t] = a[t - p];
            db(t + 1, p, k, n, a, seq);
            for j in a[t - p] + 1..k {
                a[t] = j;
                db(t + 1, t, k, n, a, seq);
            }
        }
    }
    let (mut a, mut seq) = (vec![0; k * n + 1], vec![]);
    db(1, 1, k, n, &mut a, &mut seq);
    seq
}

/// The reflected binary Gray code on n bits.
fn gray(n: u32) -> Vec<usize> {
    (0..1usize << n).map(|i| i ^ (i >> 1)).collect()
}

fn windows<T: Clone>(seq: &[T], n: usize, cyclic: bool) -> Vec<Vec<T>> {
    let s: Vec<T> = if cyclic { seq.iter().chain(&seq[..n - 1]).cloned().collect() } else { seq.to_vec() };
    let m = if cyclic { seq.len() } else { seq.len() + 1 - n };
    (0..m).map(|i| s[i..i + n].to_vec()).collect()
}

fn bitstr(v: usize, n: usize) -> String {
    format!("{v:0n$b}")
}

// ------------------------------------------------- 1. layer memory: A B C

/// The atomic-layer-memory code: 0 = A, 1 = C, a B between every bit.
fn spacer_stack(bits: &str) -> String {
    bits.chars().map(|b| if b == '1' { "C" } else { "A" }).collect::<Vec<_>>().join("B")
}

fn moved(s: &str, t: &str) -> usize {
    s.chars().zip(t.chars()).filter(|(x, y)| x != y).count()
}

/// Stacks for counting 0..2^n-1 in the given code order, and how many layers
/// move at each step.
fn counting(n: usize, codes: &[usize]) -> (Vec<String>, Vec<usize>) {
    let st: Vec<String> = codes.iter().map(|&c| spacer_stack(&bitstr(c, n))).collect();
    let mv = (0..st.len() - 1).map(|i| moved(&st[i], &st[i + 1])).collect();
    (st, mv)
}

/// The spacer code wrapped round a binary de Bruijn sequence: the B spacers
/// are ignored, and every n-bit word turns up exactly once among the bit
/// layers, read in windows of n and wrapping round.
fn spacer_de_bruijn(n: usize) -> (Vec<char>, String) {
    let bits: String = de_bruijn(2, n).iter().map(|d| d.to_string()).collect();
    ((spacer_stack(&bits) + "B").chars().collect(), bits)
}

// ------------------------------------------------- 2. the Logic Alphabet

const ROWS: [&str; 4] = ["TT", "TF", "FT", "FF"];
const MEANING: [(&str, &str); 16] = [
    ("1111", "always true"),
    ("0000", "never true"),
    ("1000", "A and B"),
    ("0111", "not both A and B"),
    ("1110", "A or B"),
    ("0001", "neither A nor B"),
    ("1011", "A implies B"),
    ("0100", "A and not B"),
    ("1101", "B implies A"),
    ("0010", "B and not A"),
    ("1100", "A"),
    ("0011", "not A"),
    ("1010", "B"),
    ("0101", "not B"),
    ("1001", "A if and only if B"),
    ("0110", "A or B but not both"),
];

fn meaning(bits: &str) -> &'static str {
    MEANING.iter().find(|(b, _)| *b == bits).unwrap().1
}

/// Every Gray cycle through the 16 truth tables (a Hamiltonian cycle on the
/// 4-cube), as directed cycles from 0000.
fn logic_cycles() -> Vec<Vec<usize>> {
    fn dfs(path: &mut Vec<usize>, seen: u32, res: &mut Vec<Vec<usize>>) {
        let v = *path.last().unwrap();
        if path.len() == 16 {
            if v.count_ones() == 1 {
                res.push(path.clone());
            }
            return;
        }
        for b in 0..4 {
            let w = v ^ (1 << b);
            if seen >> w & 1 == 0 {
                path.push(w);
                dfs(path, seen | 1 << w, res);
                path.pop();
            }
        }
    }
    let mut res = vec![];
    dfs(&mut vec![0], 1, &mut res);
    res
}

fn antipodal(c: &[usize]) -> bool {
    (0..16).all(|i| c[(i + 8) % 16] == c[i] ^ 15)
}

/// The 48 antipodal Gray cycles, each turned to start at 1111 and read the
/// way round whose second glyph is the larger number, highest first.
fn logic_rings() -> Vec<Vec<String>> {
    let mut seen = BTreeSet::new();
    for c in logic_cycles() {
        if antipodal(&c) {
            let i = c.iter().position(|&v| v == 15).unwrap();
            let mut r: Vec<usize> = c[i..].iter().chain(&c[..i]).cloned().collect();
            if r[1] < r[15] {
                r[1..].reverse();
            }
            seen.insert(r);
        }
    }
    seen.into_iter().rev().map(|r| r.iter().map(|&v| bitstr(v, 4)).collect()).collect()
}

/// The ring on the cut file: the first of the 48.
fn logic_ring() -> Vec<String> {
    logic_rings().swap_remove(0)
}

// ------------------------------------------------- 3. the grid script

const DIRS: [&str; 8] = ["E", "NE", "N", "NW", "W", "SW", "S", "SE"]; // 45 deg apart

fn step(d: &str) -> (i32, i32) {
    match d {
        "E" => (1, 0),
        "W" => (-1, 0),
        "N" => (0, -1),
        "S" => (0, 1),
        "NE" => (1, -1),
        "NW" => (-1, -1),
        "SE" => (1, 1),
        "SW" => (-1, 1),
        _ => unreachable!(),
    }
}

type Trace = (&'static str, &'static str);

/// the 28 letters as gernreich-writing/gernreich.js traces them
const LETTERS: [(&str, Trace); 28] = [
    ("KH", ("N", "SE")),
    ("Ř", ("E", "SW")),
    ("Z", ("SE", "N")),
    ("Q", ("SW", "E")),
    ("J", ("S", "NE")),
    ("X", ("W", "SE")),
    ("K", ("NE", "S")),
    ("V", ("SE", "W")),
    ("B", ("S", "E")),
    ("W", ("N", "E")),
    ("P", ("E", "S")),
    ("Y", ("S", "W")),
    ("G", ("E", "E")),
    ("F", ("S", "S")),
    ("M", ("SE", "SW")),
    ("C", ("SE", "NE")),
    ("U", ("NE", "SE")),
    ("L", ("SW", "SE")),
    ("D", ("SW", "S")),
    ("H", ("E", "SE")),
    ("S", ("S", "SW")),
    ("I", ("SE", "E")),
    ("R", ("SE", "S")),
    ("N", ("E", "NE")),
    ("O", ("S", "SE")),
    ("A", ("NE", "E")),
    ("T", ("SE", "SE")),
    ("E", ("NE", "NE")),
];

fn di(d: &str) -> i32 {
    DIRS.iter().position(|x| *x == d).unwrap() as i32
}

fn opp(d: &str) -> &'static str {
    DIRS[((di(d) + 4) % 8) as usize]
}

/// How far b is turned from a, in eighths of a turn, 0..4.
fn eighths(a: &str, b: &str) -> i32 {
    let t = (di(b) - di(a)).rem_euclid(8);
    t.min(8 - t)
}

/// Signed turn at the letter's middle, in eighths, -3..3.
fn bend(a: &str, b: &str) -> i32 {
    (di(b) - di(a) + 4).rem_euclid(8) - 4
}

/// a letter reads either way
fn tracings(p: Trace) -> [Trace; 2] {
    [p, (opp(p.1), opp(p.0))]
}

fn letter(k: usize) -> &'static str {
    LETTERS[k].0
}

fn index_of(name: &str) -> usize {
    LETTERS.iter().position(|(k, _)| *k == name).unwrap()
}

fn shape() -> HashMap<Trace, usize> {
    let mut m = HashMap::new();
    for (k, (_, p)) in LETTERS.iter().enumerate() {
        for t in tracings(*p) {
            m.insert(t, k);
        }
    }
    m
}

/// Neighbours of each letter as (letter, cost), in the order Python's dict
/// first met them.
type Moves = Vec<Vec<(usize, i32)>>;

/// Letters one change apart: some tracing of each keeps one step and turns
/// the other. The cost is that turn in eighths.
fn grid_moves() -> Moves {
    let sh = shape();
    let mut nb: Moves = vec![vec![]; LETTERS.len()];
    for (k, (_, p)) in LETTERS.iter().enumerate() {
        for (a, b) in tracings(*p) {
            for c in DIRS {
                for (q, cost) in [((c, b), eighths(a, c)), ((a, c), eighths(b, c))] {
                    if let Some(&o) = sh.get(&q) {
                        if o != k {
                            match nb[k].iter_mut().find(|e| e.0 == o) {
                                Some(e) => e.1 = e.1.min(cost),
                                None => nb[k].push((o, cost.min(9))),
                            }
                        }
                    }
                }
            }
        }
    }
    nb
}

fn cost_of(nb: &Moves, v: usize, w: usize) -> Option<i32> {
    nb[v].iter().find(|e| e.0 == w).map(|e| e.1)
}

/// Turning one step by an eighth changes the bend by one, so letters of odd
/// bend and even bend alternate along such moves.
fn grid_sides() -> (Vec<&'static str>, Vec<&'static str>) {
    let mut odd: Vec<&str> = LETTERS.iter().filter(|(_, (a, b))| bend(a, b) % 2 != 0).map(|l| l.0).collect();
    odd.sort();
    let mut even: Vec<&str> = LETTERS.iter().map(|l| l.0).filter(|k| !odd.contains(k)).collect();
    even.sort();
    (even, odd)
}

/// A cycle through all 28 crosses between the sides an even number of times,
/// at most twice the smaller side; every other move stays on its side and
/// costs at least two eighths.
fn grid_bound() -> i32 {
    let (even, odd) = grid_sides();
    let cross = 2 * even.len().min(odd.len()) as i32;
    cross + 2 * (LETTERS.len() as i32 - cross)
}

/// The least-turning search from G, trying each letter's neighbours in the
/// given order; the first ring that meets grid_bound().
fn least_ring(nb: &Moves, order: &[Vec<(usize, i32)>], bound: i32) -> Option<Vec<usize>> {
    fn dfs(nb: &Moves, order: &[Vec<(usize, i32)>], bound: i32, path: &mut Vec<usize>, used: &mut [bool], cost: i32, found: &mut Option<Vec<usize>>) {
        if found.is_some() || cost + (LETTERS.len() - path.len()) as i32 > bound {
            return;
        }
        let v = *path.last().unwrap();
        if path.len() == LETTERS.len() {
            if let Some(c) = cost_of(nb, v, path[0]) {
                if cost + c == bound {
                    *found = Some(path.clone());
                }
            }
            return;
        }
        for &(w, c) in &order[v] {
            if !used[w] {
                path.push(w);
                used[w] = true;
                dfs(nb, order, bound, path, used, cost + c, found);
                used[w] = false;
                path.pop();
            }
        }
    }
    let g = index_of("G");
    let (mut path, mut used, mut found) = (vec![g], vec![false; LETTERS.len()], None);
    used[g] = true;
    dfs(nb, order, bound, &mut path, &mut used, 0, &mut found);
    found
}

/// The first Gray cycle through the 28 letters, from G, that meets
/// grid_bound(): the least turning possible.
fn grid_ring() -> Vec<&'static str> {
    let nb = grid_moves();
    let order: Vec<Vec<(usize, i32)>> = nb
        .iter()
        .map(|v| {
            let mut v = v.clone();
            v.sort_by(|x, y| (x.1, letter(x.0)).cmp(&(y.1, letter(y.0))));
            v
        })
        .collect();
    least_ring(&nb, &order, grid_bound()).unwrap().into_iter().map(letter).collect()
}

/// A JSON list of lists of strings, as json.load reads grid-rings.json.
fn parse_string_lists(s: &str) -> Vec<Vec<String>> {
    let (mut out, mut cur, mut depth) = (vec![], vec![], 0);
    let mut it = s.chars();
    while let Some(c) = it.next() {
        match c {
            '[' => {
                depth += 1;
                if depth == 2 {
                    cur = vec![];
                }
            }
            ']' => {
                if depth == 2 {
                    out.push(std::mem::take(&mut cur));
                }
                depth -= 1;
            }
            '"' => {
                let mut t = String::new();
                loop {
                    match it.next().expect("unterminated string") {
                        '"' => break,
                        '\\' => match it.next().unwrap() {
                            'u' => {
                                let h: String = (0..4).map(|_| it.next().unwrap()).collect();
                                t.push(char::from_u32(u32::from_str_radix(&h, 16).unwrap()).unwrap());
                            }
                            e => t.push(e),
                        },
                        ch => t.push(ch),
                    }
                }
                cur.push(t);
            }
            _ => {}
        }
    }
    out
}

/// The sample of least-turning grid rings the page shows: read from
/// grid-rings.json (written by --resample).
fn grid_rings() -> Vec<Vec<String>> {
    parse_string_lists(&std::fs::read_to_string(here().join("grid-rings.json")).unwrap())
}

/// Python's random.Random(seed): MT19937 seeded by init_by_array, and
/// random() from two draws, 53 bits.
struct PyRandom {
    mt: [u32; 624],
    i: usize,
}

impl PyRandom {
    fn new(seed: u32) -> PyRandom {
        let mut mt = [0u32; 624];
        mt[0] = 19650218;
        for i in 1..624 {
            mt[i] = 1812433253u32.wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_add(i as u32);
        }
        let key = [seed];
        let (mut i, mut j) = (1usize, 0usize);
        for _ in 0..624.max(key.len()) {
            mt[i] = (mt[i] ^ (mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1664525)).wrapping_add(key[j]).wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= 624 {
                mt[0] = mt[623];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
        }
        for _ in 0..623 {
            mt[i] = (mt[i] ^ (mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1566083941)).wrapping_sub(i as u32);
            i += 1;
            if i >= 624 {
                mt[0] = mt[623];
                i = 1;
            }
        }
        mt[0] = 0x8000_0000;
        PyRandom { mt, i: 624 }
    }

    fn next_u32(&mut self) -> u32 {
        if self.i >= 624 {
            for k in 0..624 {
                let y = (self.mt[k] & 0x8000_0000) | (self.mt[(k + 1) % 624] & 0x7fff_ffff);
                self.mt[k] = self.mt[(k + 397) % 624] ^ (y >> 1) ^ if y & 1 == 1 { 0x9908_b0df } else { 0 };
            }
            self.i = 0;
        }
        let mut y = self.mt[self.i];
        self.i += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^ (y >> 18)
    }

    fn random(&mut self) -> f64 {
        let (a, b) = (self.next_u32() >> 5, self.next_u32() >> 6);
        (a as f64 * 67108864.0 + b as f64) * (1.0 / 9007199254740992.0)
    }
}

/// Several different least-turning rings from G: the first is grid_ring();
/// the rest come from the same search with the neighbours tried in a
/// shuffled order (seeded, so the list is the same every run). A ring and
/// its reverse count as one.
fn sample_grid_rings(count: usize, seed: u32) -> Vec<Vec<&'static str>> {
    let (nb, bound) = (grid_moves(), grid_bound());
    let first: Vec<usize> = grid_ring().into_iter().map(index_of).collect();
    let reverse = |r: &Vec<usize>| -> Vec<usize> { std::iter::once(r[0]).chain(r[1..].iter().rev().cloned()).collect() };
    let mut seen: HashSet<Vec<usize>> = [first.clone(), reverse(&first)].into_iter().collect();
    let mut out = vec![first];
    let mut rnd = PyRandom::new(seed);
    let mut tries = 0;
    while out.len() < count && tries < count * 20 {
        tries += 1;
        let order: Vec<Vec<(usize, i32)>> = nb
            .iter()
            .map(|v| {
                let mut keyed: Vec<((i32, f64), (usize, i32))> = v.iter().map(|&e| ((e.1, rnd.random()), e)).collect();
                keyed.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
                keyed.into_iter().map(|x| x.1).collect()
            })
            .collect();
        if let Some(r) = least_ring(&nb, &order, bound) {
            if !seen.contains(&r) {
                seen.insert(reverse(&r));
                seen.insert(r.clone());
                out.push(r);
            }
        }
    }
    out.into_iter().map(|r| r.into_iter().map(letter).collect()).collect()
}

/// How letter k becomes o: the tracing of k, the tracing of o, the slot
/// turned (0 or 1) and the turn in eighths.
fn grid_change(k: &str, o: &str) -> (Trace, Trace, usize, i32) {
    let mut best: Option<(Trace, Trace, usize, i32)> = None;
    for ab in tracings(LETTERS[index_of(k)].1) {
        for q in tracings(LETTERS[index_of(o)].1) {
            for slot in 0..2 {
                let pick = |t: Trace, s: usize| if s == 0 { t.0 } else { t.1 };
                if pick(ab, 1 - slot) == pick(q, 1 - slot) {
                    let t = eighths(pick(ab, slot), pick(q, slot));
                    if best.is_none_or(|b| t < b.3) {
                        best = Some((ab, q, slot, t));
                    }
                }
            }
        }
    }
    best.unwrap()
}

// ---------------------------------------------------------------- checks

fn check() -> bool {
    let mut out: Vec<String> = vec![];
    let mut ok = |cond: bool, what: String| out.push(format!("{}{what}", if cond { "PASS " } else { "FAIL " }));

    for (k, n) in [(2usize, 5usize), (2, 9)] {
        let s = de_bruijn(k, n);
        let w: HashSet<Vec<usize>> = windows(&s, n, true).into_iter().collect();
        let kn = k.pow(n as u32);
        ok(s.len() == kn && w.len() == kn, format!("B({k},{n}) has length {kn} and every window of {n} once"));
    }
    let track = ruler_track(300);
    ok(windows(&track, 9, false).into_iter().collect::<HashSet<_>>().len() == 301, "the ruler track: 301 windows of 9 cells, all different".into());
    for n in [5u32] {
        let g = gray(n);
        ok((0..g.len()).all(|i| (g[i] ^ g[(i + 1) % g.len()]).count_ones() == 1), format!("the {n}-bit Gray code changes one bit per step, round the circle"));
    }
    let (st, mv) = counting(4, &(0..16).collect::<Vec<_>>());
    let (gst, gmv) = counting(4, &gray(4));
    let (smv, sgmv) = (mv.iter().sum::<usize>(), gmv.iter().sum::<usize>());
    ok(smv == 26 && sgmv == 15 && *gmv.iter().max().unwrap() == 1, format!("counting 0..15 in the layer stack moves {smv} layers in binary, {sgmv} in Gray"));
    ok(
        st[0] == gst[0] && st.iter().chain(&gst).all(|s| s.as_bytes().iter().skip(1).step_by(2).all(|&c| c == b'B')),
        "the B spacers never move: only the A and C bit layers count".into(),
    );
    for n in [3usize, 4, 5] {
        let (s, _) = spacer_de_bruijn(n);
        let bitl: Vec<char> = s.iter().step_by(2).cloned().collect();
        let w: HashSet<Vec<char>> = windows(&bitl, n, true).into_iter().collect();
        let l = s.len();
        ok(
            l == 1 << (n + 1) && (0..l).all(|i| s[i] != s[(i + 1) % l]) && s.iter().skip(1).step_by(2).all(|&c| c == 'B') && s.iter().skip(1).step_by(2).count() == 1 << n && w.len() == 1 << n,
            format!("spacer stack, windows of {n} bits: {l} layers, B ignored, every {n}-bit word once"),
        );
    }
    let cyc = logic_cycles();
    ok(cyc.len() == 2688, format!("{} Gray cycles through the 16 truth tables", cyc.len() / 2));
    let anti = cyc.iter().filter(|c| antipodal(c)).count();
    ok(anti == 96, format!("{} of them antipodal", anti / 2));
    let mut all16: Vec<&str> = MEANING.iter().map(|m| m.0).collect();
    all16.sort();
    fn sorted(r: &[String]) -> Vec<&str> {
        let mut v: Vec<&str> = r.iter().map(|s| s.as_str()).collect();
        v.sort();
        v
    }
    let one_step = |r: &Vec<String>| (0..16).all(|i| moved(&r[i], &r[(i + 1) % 16]) == 1);
    let rings = logic_rings();
    ok(
        rings.len() == 48
            && rings.iter().collect::<HashSet<_>>().len() == 48
            && rings.iter().all(|r| sorted(r) == all16 && r[0] == "1111" && antipodal(&r.iter().map(|t| usize::from_str_radix(t, 2).unwrap()).collect::<Vec<_>>()) && one_step(r)),
        "the 48 rings listed: all different, all antipodal Gray cycles from 1111".into(),
    );
    let ring = logic_ring();
    let num = |t: &str| usize::from_str_radix(t, 2).unwrap();
    ok(
        sorted(&ring) == all16 && one_step(&ring) && (0..8).all(|i| num(&ring[i]) ^ num(&ring[i + 8]) == 15),
        "the Logic ring: all 16, one square per step, each across from its negation".into(),
    );
    let sh = shape();
    ok(sh.len() == 56 && sh.values().collect::<HashSet<_>>().len() == 28, "28 letters, 56 tracings".into());
    let nb = grid_moves();
    let mut all28: Vec<&str> = LETTERS.iter().map(|l| l.0).collect();
    all28.sort();
    let grs = grid_rings();
    ok(
        grs.len() == 48
            && grs.iter().collect::<HashSet<_>>().len() == 48
            && grs.iter().all(|r| {
                let mut s: Vec<&str> = r.iter().map(|x| x.as_str()).collect();
                s.sort();
                s == all28 && r[0] == "G" && (0..28).map(|i| cost_of(&nb, index_of(&r[i]), index_of(&r[(i + 1) % 28])).unwrap()).sum::<i32>() == 32
            }),
        "48 different least-turning grid rings listed, each 32 eighths".into(),
    );
    let (even, odd) = grid_sides();
    ok((even.len(), odd.len()) == (12, 16), format!("grid letters: {} of even bend, {} of odd", even.len(), odd.len()));
    ok(
        nb.iter().enumerate().all(|(k, v)| v.iter().filter(|e| e.1 == 1).all(|e| even.contains(&letter(k)) != even.contains(&letter(e.0)))),
        "every one-eighth move crosses between the sides, so no cycle uses only those".into(),
    );
    let gr = grid_ring();
    let cost: Vec<i32> = (0..28).map(|i| cost_of(&nb, index_of(gr[i]), index_of(gr[(i + 1) % 28])).unwrap()).collect();
    let mut gs = gr.clone();
    gs.sort();
    let total: i32 = cost.iter().sum();
    ok(
        gs == all28 && total == grid_bound() && grid_bound() == 32,
        format!(
            "the grid ring: all 28, {} one-eighth turns and {} quarter turns, {total} eighths, the least possible",
            cost.iter().filter(|&&c| c == 1).count(),
            cost.iter().filter(|&&c| c == 2).count()
        ),
    );
    for line in &out {
        println!("{line}");
    }
    out.iter().all(|l| l.starts_with("PASS"))
}

// ------------------------------------------------------------ cut files

/// A number to 3 places, trailing zeros dropped, as walks.py's f().
fn f(v: f64) -> String {
    let s = format!("{v:.3}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    if s == "-0" { "0".into() } else { s.into() }
}

/// Python's str() of a float: shortest round trip, with ".0" on whole numbers.
fn py_float(v: f64) -> String {
    if v.fract() == 0.0 && v.abs() < 1e16 { format!("{v:.1}") } else { format!("{v}") }
}

type Pt = (f64, f64);

/// single-stroke digits on a 1 x 2 box, y down (from staircase_ruler.py)
const DIGITS: [(char, &[&[Pt]]); 10] = [
    ('0', &[&[(0.2, 0.0), (0.8, 0.0), (1.0, 0.3), (1.0, 1.7), (0.8, 2.0), (0.2, 2.0), (0.0, 1.7), (0.0, 0.3), (0.2, 0.0)]]),
    ('1', &[&[(0.25, 0.4), (0.6, 0.0), (0.6, 2.0)]]),
    ('2', &[&[(0.0, 0.4), (0.25, 0.0), (0.75, 0.0), (1.0, 0.4), (1.0, 0.8), (0.0, 2.0), (1.0, 2.0)]]),
    ('3', &[&[(0.0, 0.0), (1.0, 0.0), (0.4, 0.8), (0.75, 0.8), (1.0, 1.1), (1.0, 1.7), (0.75, 2.0), (0.25, 2.0), (0.0, 1.7)]]),
    ('4', &[&[(0.75, 2.0), (0.75, 0.0), (0.0, 1.4), (1.0, 1.4)]]),
    ('5', &[&[(1.0, 0.0), (0.1, 0.0), (0.0, 0.9), (0.7, 0.8), (1.0, 1.1), (1.0, 1.7), (0.75, 2.0), (0.25, 2.0), (0.0, 1.7)]]),
    ('6', &[&[(0.9, 0.0), (0.4, 0.0), (0.0, 0.6), (0.0, 1.7), (0.25, 2.0), (0.75, 2.0), (1.0, 1.7), (1.0, 1.2), (0.75, 0.9), (0.25, 0.9), (0.0, 1.2)]]),
    ('7', &[&[(0.0, 0.0), (1.0, 0.0), (0.35, 2.0)]]),
    ('8', &[&[(0.5, 0.9), (0.1, 0.7), (0.1, 0.25), (0.35, 0.0), (0.65, 0.0), (0.9, 0.25), (0.9, 0.7), (0.5, 0.9), (0.0, 1.2), (0.0, 1.7), (0.25, 2.0), (0.75, 2.0), (1.0, 1.7), (1.0, 1.2), (0.5, 0.9)]]),
    ('9', &[&[(0.1, 2.0), (0.6, 2.0), (1.0, 1.4), (1.0, 0.3), (0.75, 0.0), (0.25, 0.0), (0.0, 0.3), (0.0, 0.8), (0.25, 1.1), (0.75, 1.1), (1.0, 0.8)]]),
];

fn number(text: &str, cx: f64, cy: f64, h: f64) -> Vec<Vec<Pt>> {
    let (w, gap) = (h / 2.0, h * 0.15);
    let len = text.chars().count() as f64;
    let x0 = cx - (len * w + (len - 1.0) * gap) / 2.0;
    let mut out = vec![];
    for (i, ch) in text.chars().enumerate() {
        for s in DIGITS.iter().find(|d| d.0 == ch).unwrap().1.iter() {
            out.push(s.iter().map(|&(px, py)| (x0 + i as f64 * (w + gap) + px * w, cy - h / 2.0 + py * w)).collect());
        }
    }
    out
}

fn pl(pts: &[Pt], close: bool) -> String {
    "M".to_string() + &pts.iter().map(|&(x, y)| format!("{},{}", f(x), f(y))).collect::<Vec<_>>().join(" L") + if close { " Z" } else { "" }
}

fn circle(cx: f64, cy: f64, r: f64) -> String {
    format!("M{},{} A{},{} 0 1 0 {},{} A{},{} 0 1 0 {},{} Z", f(cx - r), f(cy), f(r), f(r), f(cx + r), f(cy), f(r), f(r), f(cx - r), f(cy))
}

fn rect(x: f64, y: f64, w: f64, h: f64) -> String {
    format!("M{},{} H{} V{} H{} Z", f(x), f(y), f(x + w), f(y + h), f(x))
}

/// Point at radius r, `turn` of a full turn clockwise from the top.
fn polar(cx: f64, cy: f64, r: f64, turn: f64) -> Pt {
    let a = 2.0 * PI * turn;
    (cx + r * a.sin(), cy - r * a.cos())
}

/// Annular sector from turn t0 to t1 clockwise, radii r0 < r1.
fn wedge(cx: f64, cy: f64, r0: f64, r1: f64, t0: f64, t1: f64) -> String {
    let big = if t1 - t0 > 0.5 { 1 } else { 0 };
    let (a, b) = (polar(cx, cy, r1, t0), polar(cx, cy, r1, t1));
    let (c, d) = (polar(cx, cy, r0, t1), polar(cx, cy, r0, t0));
    format!(
        "M{},{} A{},{} 0 {big} 1 {},{} L{},{} A{},{} 0 {big} 0 {},{} Z",
        f(a.0),
        f(a.1),
        f(r1),
        f(r1),
        f(b.0),
        f(b.1),
        f(c.0),
        f(c.1),
        f(r0),
        f(r0),
        f(d.0),
        f(d.1)
    )
}

/// (start, length) of each run of 1s in a list.
fn runs(cells: &[usize]) -> Vec<(usize, usize)> {
    let (mut res, mut i) = (vec![], 0);
    while i < cells.len() {
        if cells[i] != 0 {
            let mut j = i;
            while j < cells.len() && cells[j] != 0 {
                j += 1;
            }
            res.push((i, j - i));
            i = j;
        } else {
            i += 1;
        }
    }
    res
}

/// Paths by layer; written as a cut file or as a page preview.
struct Piece {
    name: String,
    bx: (f64, f64, f64, f64), // x0, y0, x1, y1
    engrave: Vec<String>,
    score: Vec<String>,
    inner: Vec<String>,
    outline: Vec<String>,
}

impl Piece {
    fn new(name: &str, bx: (f64, f64, f64, f64)) -> Piece {
        Piece { name: name.into(), bx, engrave: vec![], score: vec![], inner: vec![], outline: vec![] }
    }

    fn size(&self) -> (f64, f64) {
        (self.bx.2 - self.bx.0, self.bx.3 - self.bx.1)
    }

    fn cut_svg(&self) -> String {
        let (x0, y0, x1, y1) = self.bx;
        let p = 1.0;
        let (w, h) = (x1 - x0 + 2.0 * p, y1 - y0 + 2.0 * p);
        let paths = |ds: &[String]| ds.iter().map(|d| format!("    <path d=\"{d}\"/>\n")).collect::<String>();
        let g = |col: &str, ds: &[String]| format!("  <g fill=\"none\" stroke=\"{col}\" stroke-width=\"0.1\">\n{}  </g>\n", paths(ds));
        format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{}mm\" height=\"{}mm\" viewBox=\"{} {} {} {}\">\n  <g fill=\"{GREEN}\" stroke=\"none\">\n{}  </g>\n{}{}{}</svg>\n",
            f(w),
            f(h),
            f(x0 - p),
            f(y0 - p),
            f(w),
            f(h),
            paths(&self.engrave),
            g(BLUE, &self.score),
            g(ORANGE, &self.inner),
            g(BLACK, &self.outline)
        )
    }

    fn preview(&self, label: &str) -> String {
        let (x0, y0, x1, y1) = self.bx;
        let p = 1.0;
        let (w, h) = (x1 - x0 + 2.0 * p, y1 - y0 + 2.0 * p);
        let body = self.outline.iter().chain(&self.inner).cloned().collect::<Vec<_>>().join(" ");
        format!(
            "<svg class=\"piece\" viewBox=\"{} {} {} {}\" role=\"img\" aria-label=\"{label}\"><path class=\"stock\" fill-rule=\"evenodd\" d=\"{body}\"/>{}<path class=\"sco\" d=\"{}\"/></svg>",
            f(x0 - p),
            f(y0 - p),
            f(w),
            f(h),
            self.engrave.iter().map(|d| format!("<path class=\"eng\" d=\"{d}\"/>")).collect::<String>(),
            self.score.join(" ")
        )
    }
}

// 1. the window ruler -------------------------------------------------------

const RULER_W: f64 = 30.0;
const TRACK_H: f64 = 8.0;
const TAIL: f64 = 20.0;
const HOLE_D: f64 = 6.0;
const CORNER_R: f64 = 4.0;
const BLEED: f64 = 0.3;
const WINDOW: usize = 9;

/// Cells of the ruler's track, 1mm each: the start of B(2,9) unrolled, long
/// enough that a 9mm window at every whole mm 0..length sits on it.
fn ruler_track(length: usize) -> Vec<usize> {
    let s = de_bruijn(2, WINDOW);
    s.iter().chain(&s[..WINDOW - 1]).take(length + WINDOW).cloned().collect()
}

fn window_ruler(length: usize, kerf: f64) -> Piece {
    let k = kerf / 2.0;
    let w = RULER_W;
    let lf = length as f64;
    let mut pc = Piece::new(&format!("window-ruler-{length}"), (-k, -k, lf + TAIL + k, w + k));
    let track = ruler_track(length);
    let ty = w - TRACK_H;
    for (x, n) in runs(&track) {
        pc.engrave.push(rect(x as f64, ty, n as f64, TRACK_H + BLEED));
    }
    for mm in 0..=track.len() {
        let m = mm as f64;
        pc.score.push(pl(&[(m, ty), (m, ty - if mm % 10 == 0 { 3.0 } else { 1.5 })], false));
    }
    for mm in 1..=length {
        let n = mm % 10;
        let m = mm as f64;
        pc.score.push(pl(&[(m, 0.0), (m, if n == 0 { 7.0 } else if n == 5 { 5.0 } else { 3.0 })], false));
    }
    for mm in (10..=length).step_by(10) {
        pc.score.extend(number(&mm.to_string(), mm as f64, 11.5, 4.0).iter().map(|s| pl(s, false)));
    }
    let (x0, y0, x1, y1, r) = (-k, -k, lf + TAIL + k, w + k, CORNER_R + k);
    pc.outline.push(format!(
        "M{},{} H{} A{},{} 0 0 1 {},{} V{} A{},{} 0 0 1 {},{} H{} Z",
        f(x0),
        f(y0),
        f(x1 - r),
        f(r),
        f(r),
        f(x1),
        f(y0 + r),
        f(y1 - r),
        f(r),
        f(r),
        f(x1 - r),
        f(y1),
        f(x0)
    ));
    pc.inner.push(circle(lf + TAIL / 2.0, w / 2.0, HOLE_D / 2.0 - k));
    pc
}

/// A card to lay on the track: a window 9 cells wide, ticks at every cell
/// edge, and a line from the window's left edge to read the scale by.
fn window_reader(kerf: f64) -> Piece {
    let k = kerf / 2.0;
    let (w, h) = (36.0, 30.0);
    let win = WINDOW as f64;
    let mut pc = Piece::new("window-reader", (-k, -k, w + k, h + k));
    let (wx, wy) = ((w - win) / 2.0, h - 4.0 - TRACK_H);
    pc.inner.push(rect(wx + k, wy + k, win - 2.0 * k, TRACK_H - 2.0 * k));
    for i in 0..=WINDOW {
        let x = wx + i as f64;
        pc.score.push(pl(&[(x, wy - 0.6), (x, wy - if i == 0 || i == WINDOW { 2.5 } else { 1.5 })], false));
    }
    pc.score.push(pl(&[(wx, wy - 3.0), (wx, 0.0)], false));
    pc.score.push(pl(&[(wx, wy + TRACK_H + 0.6), (wx, h)], false));
    pc.outline.push(rect(-k, -k, w + 2.0 * k, h + 2.0 * k));
    pc
}

// 2. two encoder discs --------------------------------------------------------

const DISC_R: f64 = 45.0;
const DISC_HOLE: f64 = 6.0;

fn disc_ticks(pc: &mut Piece, cx: f64, cy: f64, n: usize, inner_ticks: bool) {
    for i in 0..n {
        let t = i as f64 / n as f64;
        let (a, b) = (polar(cx, cy, 42.6, t), polar(cx, cy, if i % 8 != 0 { 44.0 } else { 44.4 }, t));
        pc.score.push(pl(&[a, b], false));
        if inner_ticks {
            pc.score.push(pl(&[polar(cx, cy, 28.0, t), polar(cx, cy, 29.4, t)], false));
        }
    }
}

fn gray_disc(kerf: f64) -> Piece {
    let (k, cx, cy) = (kerf / 2.0, 0.0, 0.0);
    let mut pc = Piece::new("gray-disc", (cx - DISC_R - k, cy - DISC_R - k, cx + DISC_R + k, cy + DISC_R + k));
    let (g, n) = (gray(5), 32usize);
    for ring in 0..5 {
        let (r0, r1) = ((12 + 6 * ring) as f64, (18 + 6 * ring) as f64);
        let cells: Vec<usize> = (0..n).map(|i| (g[i] >> (4 - ring)) & 1).collect();
        for (s, m) in runs(&cells) {
            pc.engrave.push(wedge(cx, cy, r0, r1, s as f64 / n as f64, (s + m) as f64 / n as f64));
        }
        pc.score.push(if ring == 4 { circle(cx, cy, r1) } else { circle(cx, cy, r0) });
    }
    disc_ticks(&mut pc, cx, cy, n, false);
    pc.inner.push(circle(cx, cy, DISC_HOLE / 2.0 - k));
    pc.outline.push(circle(cx, cy, DISC_R + k));
    pc
}

fn de_bruijn_disc(kerf: f64) -> Piece {
    let (k, cx, cy) = (kerf / 2.0, 0.0, 0.0);
    let mut pc = Piece::new("de-bruijn-disc", (cx - DISC_R - k, cy - DISC_R - k, cx + DISC_R + k, cy + DISC_R + k));
    let (s, n) = (de_bruijn(2, 5), 32usize);
    for (st, m) in runs(&s) {
        pc.engrave.push(wedge(cx, cy, 30.0, 42.0, st as f64 / n as f64, (st + m) as f64 / n as f64));
    }
    pc.score.push(circle(cx, cy, 30.0));
    pc.score.push(circle(cx, cy, 42.0));
    disc_ticks(&mut pc, cx, cy, n, true);
    pc.inner.push(circle(cx, cy, DISC_HOLE / 2.0 - k));
    pc.outline.push(circle(cx, cy, DISC_R + k));
    pc
}

// 3. the Logic ring and 4. the grid ring --------------------------------------

const LR_R: f64 = 62.0;
const LR_SMALL: f64 = 40.0;
const LR_C: f64 = 51.0;
const CELL: f64 = 5.0;
const GAP: f64 = 1.0;

fn logic_ring_piece(kerf: f64) -> Piece {
    let k = kerf / 2.0;
    let mut pc = Piece::new("logic-ring", (-LR_R - k, -LR_R - k, LR_R + k, LR_R + k));
    for (i, t) in logic_ring().iter().enumerate() {
        let (cx, cy) = polar(0.0, 0.0, LR_C, i as f64 / 16.0);
        for (j, b) in t.chars().enumerate() {
            let x = cx - CELL - GAP / 2.0 + ((j % 2) as f64) * (CELL + GAP);
            let y = cy - CELL - GAP / 2.0 + ((j / 2) as f64) * (CELL + GAP);
            pc.score.push(rect(x, y, CELL, CELL));
            if b == '1' {
                pc.engrave.push(rect(x + 0.4, y + 0.4, CELL - 0.8, CELL - 0.8));
            }
        }
    }
    // this way up: the glyphs only mean what they mean upright
    pc.score.push(pl(&[(-2.0, -LR_SMALL - 3.5), (0.0, -LR_SMALL - 6.0), (2.0, -LR_SMALL - 3.5)], false));
    pc.inner.push(circle(0.0, 0.0, LR_SMALL - k));
    pc.outline.push(circle(0.0, 0.0, LR_R + k));
    pc
}

const GR_R: f64 = 68.0;
const GR_SMALL: f64 = 44.0;
const GR_C: f64 = 56.0;
const GSTEP: f64 = 2.6;

fn letter_points(trace: Trace, x: f64, y: f64, s: f64) -> Vec<Pt> {
    let mut pts = vec![(x, y)];
    for d in [trace.0, trace.1] {
        let (dx, dy) = step(d);
        let last = *pts.last().unwrap();
        pts.push((last.0 + s * dx as f64, last.1 + s * dy as f64));
    }
    pts
}

fn span(v: impl Iterator<Item = f64> + Clone) -> (f64, f64) {
    (v.clone().fold(f64::INFINITY, f64::min), v.fold(f64::NEG_INFINITY, f64::max))
}

fn centred(trace: Trace, cx: f64, cy: f64) -> Vec<Pt> {
    let pts = letter_points(trace, 0.0, 0.0, GSTEP);
    let (x0, x1) = span(pts.iter().map(|p| p.0));
    let (y0, y1) = span(pts.iter().map(|p| p.1));
    let (ox, oy) = (cx - (x0 + x1) / 2.0, cy - (y0 + y1) / 2.0);
    pts.iter().map(|&(px, py)| (px + ox, py + oy)).collect()
}

fn grid_ring_piece(kerf: f64) -> Piece {
    let k = kerf / 2.0;
    let mut pc = Piece::new("grid-ring", (-GR_R - k, -GR_R - k, GR_R + k, GR_R + k));
    for (i, name) in grid_ring().iter().enumerate() {
        let (cx, cy) = polar(0.0, 0.0, GR_C, i as f64 / 28.0);
        let pts = centred(LETTERS[index_of(name)].1, cx, cy);
        pc.score.push(pl(&pts, false));
        pc.engrave.push(circle(pts[0].0, pts[0].1, 0.6));
    }
    pc.score.push(pl(&[(-2.0, -GR_SMALL - 3.5), (0.0, -GR_SMALL - 6.0), (2.0, -GR_SMALL - 3.5)], false));
    pc.inner.push(circle(0.0, 0.0, GR_SMALL - k));
    pc.outline.push(circle(0.0, 0.0, GR_R + k));
    pc
}

fn pieces(kerf: f64) -> Vec<Piece> {
    vec![window_ruler(300, kerf), window_reader(kerf), gray_disc(kerf), de_bruijn_disc(kerf), logic_ring_piece(kerf), grid_ring_piece(kerf)]
}

// ---------------------------------------------------------------- the page

const HI: &str = " class=\"hi\"";
const TURNED: &str = " class=\"turned\"";

fn col(c: char) -> &'static str {
    match c {
        'A' => "var(--a)",
        'B' => "var(--b)",
        'C' => "var(--c)",
        _ => unreachable!(),
    }
}

/// One stack as a column, bottom layer at the bottom.
fn stack_svg(st: &[char], hi: &[usize], label: &str) -> String {
    let (w, h) = (14usize, 9usize);
    let n = st.len();
    let mut cells = String::new();
    for (i, &c) in st.iter().enumerate() {
        let y = 1 + (n - 1 - i) * h;
        write!(
            cells,
            "<rect x=\"1\" y=\"{y}\" width=\"{w}\" height=\"{}\" rx=\"1.5\" fill=\"{}\"{}/><text x=\"{}\" y=\"{}\">{c}</text>",
            h - 1,
            col(c),
            if hi.contains(&i) { HI } else { "" },
            py_float(1.0 + w as f64 / 2.0),
            py_float(y as f64 + h as f64 / 2.0 + 0.5)
        )
        .unwrap();
    }
    format!("<svg class=\"stack\" viewBox=\"0 0 {} {}\" width=\"{}\" height=\"{}\" role=\"img\" aria-label=\"{label}\">{cells}</svg>", w + 2, n * h + 1, w + 2, n * h + 1)
}

fn codings() -> [(&'static str, Vec<usize>); 2] {
    [("Binary", (0..16).collect()), ("Gray", gray(4))]
}

fn counting_html(n: usize) -> String {
    let mut rows = vec![];
    for (name, codes) in codings() {
        let (st, mv) = counting(n, &codes);
        let mut cols = String::new();
        for (i, s) in st.iter().enumerate() {
            let sc: Vec<char> = s.chars().collect();
            let hi: Vec<usize> = if i == 0 { vec![] } else { (0..sc.len()).filter(|&j| s.as_bytes()[j] != st[i - 1].as_bytes()[j]).collect() };
            let b = bitstr(codes[i], n);
            write!(cols, "<div class=\"col\">{}<code>{b}</code><span class=\"mv\">{}</span></div>", stack_svg(&sc, &hi, &b), if i == 0 { String::new() } else { mv[i - 1].to_string() }).unwrap();
        }
        rows.push(format!("<h3>{name}: {} layers moved</h3><div class=\"stacks\">{cols}</div>", mv.iter().sum::<usize>()));
    }
    rows.join("\n")
}

/// The counting stacks twice: as built, and with the B spacers left out.
fn counting_both(n: usize) -> String {
    let mut nob = String::new();
    for (name, codes) in codings() {
        let (st, mv) = counting(n, &codes);
        let mut cols = String::new();
        for (i, s) in st.iter().enumerate() {
            let sb = s.as_bytes();
            let hi: Vec<usize> = if i == 0 { vec![] } else { (0..sb.len()).step_by(2).filter(|&j| sb[j] != st[i - 1].as_bytes()[j]).map(|j| j / 2).collect() };
            let bits: Vec<char> = s.chars().step_by(2).collect();
            let b = bitstr(codes[i], n);
            write!(cols, "<div class=\"col\">{}<code>{b}</code><span class=\"mv\">{}</span></div>", stack_svg(&bits, &hi, &b), if i == 0 { String::new() } else { mv[i - 1].to_string() }).unwrap();
        }
        write!(nob, "<h3>{name}: {} layers moved</h3><div class=\"stacks\">{cols}</div>", mv.iter().sum::<usize>()).unwrap();
    }
    format!(
        "<p><button type=\"button\" id=\"hideB\" class=\"tog\" aria-pressed=\"false\">Hide the B layers</button></p><div class=\"withB\">{}</div><div class=\"noB\" hidden>{nob}</div>",
        counting_html(n)
    )
}

fn stack_db_html(n: usize) -> (String, Vec<char>, String) {
    let (s, bits) = spacer_de_bruijn(n);
    let sp = |c: char| if c == 'B' { " sp" } else { "" };
    let cells: String = s.iter().map(|&c| format!("<span class=\"L{}\" style=\"background:{}\">{c}</span>", sp(c), col(c))).collect();
    let wraps: String = s[..4].iter().map(|&c| format!("<span class=\"L wrap{}\">{c}</span>", sp(c))).collect();
    (format!("<div class=\"row-stack\" aria-label=\"a cyclic stack of {} layers\">{cells}{wraps}</div>", s.len()), s, bits)
}

// ------------------------------------------------- the Venn ladder

const VR: f64 = 30.0; // circle radius and centre offset of the two-set Venn
const VD: i32 = 15;
const VBOX: (i32, i32, i32, i32) = (-52, -40, 104, 80); // the universe rectangle

/// The four regions of a two-set Venn diagram as SVG paths, in truth-table
/// order: TT (A and B), TF (A only), FT (B only), FF (outside).
fn venn_regions() -> [String; 4] {
    let vh = (VR * VR - (VD * VD) as f64).sqrt();
    let (t, b, r) = (format!("0,{}", f(-vh)), format!("0,{}", f(vh)), f(VR));
    let (x, y, w, h) = VBOX;
    [
        format!("M{t} A{r},{r} 0 0 1 {b} A{r},{r} 0 0 1 {t} Z"),
        format!("M{t} A{r},{r} 0 1 0 {b} A{r},{r} 0 0 1 {t} Z"),
        format!("M{t} A{r},{r} 0 1 1 {b} A{r},{r} 0 0 0 {t} Z"),
        format!("M{x},{y} H{} V{} H{x} Z M{t} A{r},{r} 0 1 0 {b} A{r},{r} 0 1 0 {t} Z", x + w, y + h),
    ]
}

enum Labels {
    No,
    Yes,
    Rows,
}

/// A Venn diagram; each region shaded where bits has a 1, or filled with its
/// own colour from colours.
fn venn_svg(bits: &str, colours: Option<&[&str]>, labels: Labels, cls: &str) -> String {
    let (x, y, w, h) = VBOX;
    let mut out = format!("<svg class=\"{cls}\" viewBox=\"{} {} {} {}\" aria-hidden=\"true\">", x - 2, y - 2, w + 4, h + 4);
    for (i, d) in venn_regions().iter().enumerate() {
        let fill = match colours {
            Some(c) => c[i],
            None => {
                if bits.as_bytes()[i] == b'1' {
                    "var(--von)"
                } else {
                    "var(--voff)"
                }
            }
        };
        write!(out, "<path d=\"{d}\" fill=\"{fill}\" fill-rule=\"evenodd\"/>").unwrap();
    }
    write!(
        out,
        "<rect class=\"vline\" x=\"{x}\" y=\"{y}\" width=\"{w}\" height=\"{h}\"/><circle class=\"vline\" cx=\"{}\" cy=\"0\" r=\"{VR}\"/><circle class=\"vline\" cx=\"{VD}\" cy=\"0\" r=\"{VR}\"/>",
        -VD
    )
    .unwrap();
    if !matches!(labels, Labels::No) {
        out += "<text x=\"-36\" y=\"-26\">A</text><text x=\"36\" y=\"-26\">B</text>";
        if matches!(labels, Labels::Rows) {
            for (name, (lx, ly)) in ROWS.iter().zip([(0, 4), (-26, 4), (26, 4), (-44, 34)]) {
                write!(out, "<text class=\"small\" x=\"{lx}\" y=\"{ly}\">{name}</text>").unwrap();
            }
        }
    }
    out + "</svg>"
}

fn glyph_svg(bits: &str, colours: Option<&[&str]>) -> String {
    let mut cells = String::new();
    for (j, b) in bits.chars().enumerate() {
        let (x, y) = ((j % 2) * 11, (j / 2) * 11);
        let fill = match colours {
            Some(c) => c[j],
            None => {
                if b == '1' {
                    "#2b2116"
                } else {
                    "#e6d3ad"
                }
            }
        };
        write!(cells, "<rect x=\"{}\" y=\"{}\" width=\"10\" height=\"10\" rx=\"1\" fill=\"{fill}\"/>", py_float(x as f64 + 0.5), py_float(y as f64 + 0.5)).unwrap();
    }
    format!("<svg class=\"vglyph\" viewBox=\"0 0 22 22\" aria-hidden=\"true\">{cells}</svg>")
}

/// The five steps, the 16 cards and the flat 4-cube, as HTML.
fn venn_ladder() -> (String, String, String) {
    let rc = ["var(--r0)", "var(--r1)", "var(--r2)", "var(--r3)"];
    let steps: [(&str, String, &str); 5] = [
        ("Two sets", venn_svg("0000", None, Labels::Yes, "venn"), "A circle for A, a circle for B, inside a box for everything else."),
        (
            "Four regions",
            venn_svg("", Some(&rc), Labels::Rows, "venn") + &glyph_svg("1111", Some(&rc)),
            "They cut the box into four regions, one per row of the truth table, and each matches a place in the glyph square.",
        ),
        (
            "2<sup>4</sup> = 16 shadings",
            format!("<div class=\"mini16\">{}</div>", (0..16).rev().map(|v| venn_svg(&bitstr(v, 4), None, Labels::No, "venn tiny")).collect::<String>()),
            "Shade each region or leave it, independently: 2 × 2 × 2 × 2 = 16 ways.",
        ),
        (
            "16 functions",
            format!("<div class=\"fnlist\">{}<span>…</span></div>", [8, 14, 11, 9].iter().map(|&v| format!("<span>{}</span>", meaning(&bitstr(v, 4)))).collect::<String>()),
            "Each shading is a function of A and B: true in the shaded regions, false elsewhere. There are exactly 16.",
        ),
        (
            "16 corners of a 4-cube",
            tesseract_svg(true),
            "Give each region its own axis, 0 for unshaded and 1 for shaded. The 16 shadings are then the corners of a four-dimensional cube, and each edge changes one region.",
        ),
    ];
    let ladder: String = steps.iter().map(|(t, fig, txt)| format!("<li><b>{t}</b>{fig}<span>{txt}</span></li>")).collect();
    let cards: String = (0..16)
        .rev()
        .map(|v| {
            let b = bitstr(v, 4);
            format!(
                "<button type=\"button\" class=\"vcard\" data-v=\"{v}\" aria-pressed=\"false\">{}{}<span>{}</span></button>",
                venn_svg(&b, None, Labels::No, "venn"),
                glyph_svg(&b, None),
                meaning(&b)
            )
        })
        .collect();
    (format!("<ol class=\"ladder\">{ladder}</ol>"), cards, tesseract_svg(false))
}

/// The 4-cube drawn flat as a square inside a square: TT and TF are the big
/// square's sides, FT and FF a half-size square turned 25 degrees (the turn
/// and size that keep corners and edges furthest apart).
fn tesseract_svg(small: bool) -> String {
    let s = 62.0f64;
    let mut vec: Vec<Pt> = vec![(s, 0.0), (0.0, s)];
    for a in [25.0f64, 115.0] {
        vec.push((s / 2.0 * a.to_radians().cos(), s / 2.0 * a.to_radians().sin()));
    }
    let at = |v: usize| -> Pt {
        let on: Vec<usize> = (0..4).filter(|i| v >> (3 - i) & 1 == 1).collect();
        let x = on.iter().fold(0.0, |t, &i| t + vec[i].0);
        let y = -on.iter().fold(0.0, |t, &i| t + vec[i].1);
        (x, y)
    };
    let raw: Vec<Pt> = (0..16).map(at).collect();
    let (x0, x1) = span(raw.iter().map(|p| p.0));
    let (y0, y1) = span(raw.iter().map(|p| p.1));
    let (cx, cy) = ((x0 + x1) / 2.0, (y0 + y1) / 2.0);
    let pts: Vec<Pt> = raw.iter().map(|&(x, y)| (x - cx, y - cy)).collect();
    let m = pts.iter().map(|p| p.0.abs().max(p.1.abs())).fold(f64::NEG_INFINITY, f64::max) + 11.0;
    let vb = format!("{} {} {} {}", f(-m), f(-m), f(2.0 * m), f(2.0 * m));
    let mut edges = String::new();
    for v in 0..16usize {
        for b in 0..4 {
            let u = v ^ (1 << b);
            if v < u {
                write!(edges, "<line class=\"tedge\" data-a=\"{v}\" data-b=\"{u}\" x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\"/>", f(pts[v].0), f(pts[v].1), f(pts[u].0), f(pts[u].1)).unwrap();
            }
        }
    }
    if small {
        let dots: String = pts.iter().map(|&(x, y)| format!("<circle cx=\"{}\" cy=\"{}\" r=\"4\"/>", f(x), f(y))).collect();
        return format!("<svg class=\"tess small\" viewBox=\"{vb}\" aria-hidden=\"true\">{edges}{dots}</svg>");
    }
    let mut nodes = String::new();
    for (v, &(x, y)) in pts.iter().enumerate() {
        write!(nodes, "<g class=\"tnode\" data-v=\"{v}\" transform=\"translate({},{}) scale(0.7273)\"><rect x=\"-1\" y=\"-1\" width=\"24\" height=\"24\" rx=\"3\" class=\"tbg\"/>", f(x - 8.0), f(y - 8.0)).unwrap();
        for (j, b) in bitstr(v, 4).chars().enumerate() {
            write!(
                nodes,
                "<rect x=\"{}\" y=\"{}\" width=\"10\" height=\"10\" rx=\"1\" class=\"{}\"/>",
                py_float(((j % 2) * 11) as f64 + 0.5),
                py_float(((j / 2) * 11) as f64 + 0.5),
                if b == '1' { "on" } else { "off" }
            )
            .unwrap();
        }
        nodes += "</g>";
    }
    format!("<svg class=\"tess\" viewBox=\"{vb}\" role=\"img\" aria-label=\"the 16 functions as the corners of a 4-cube drawn flat\">{edges}{nodes}</svg>")
}

fn letter_svg(trace: Trace, hi_slot: usize) -> String {
    let pts = letter_points(trace, 0.0, 0.0, 10.0);
    let (xa, xb) = span(pts.iter().map(|p| p.0));
    let (ya, yb) = span(pts.iter().map(|p| p.1));
    let (x0, y0) = (xa - 4.0, ya - 4.0);
    let (w, h) = (xb - xa + 8.0, yb - ya + 8.0);
    let seg: String = (0..2)
        .map(|i| format!("<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\"{}/>", f(pts[i].0), f(pts[i].1), f(pts[i + 1].0), f(pts[i + 1].1), if i == hi_slot { TURNED } else { "" }))
        .collect();
    format!(
        "<svg class=\"letter\" viewBox=\"{} {} {} {}\" width=\"{}\" height=\"{}\" aria-hidden=\"true\">{seg}<circle cx=\"{}\" cy=\"{}\" r=\"2\"/></svg>",
        f(x0),
        f(y0),
        f(w),
        f(h),
        f(w * 2.2),
        f(h * 2.2),
        f(pts[0].0),
        f(pts[0].1)
    )
}

fn grid_list(ring: &[&str]) -> String {
    let mut out = String::new();
    for (i, name) in ring.iter().enumerate() {
        let prev = ring[(i + ring.len() - 1) % ring.len()];
        let (_, q, slot, t) = grid_change(prev, name);
        write!(out, "<div class=\"glet\">{}<b>{name}</b><span>{}</span></div>", letter_svg(q, slot), if t == 1 { "⅛" } else { "¼" }).unwrap();
    }
    out
}

/// A JSON string as Python's json.dumps writes it.
fn js(s: &str, ascii: bool) -> String {
    let mut o = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => o += "\\\"",
            '\\' => o += "\\\\",
            '\n' => o += "\\n",
            c if (c as u32) < 0x20 || (ascii && (c as u32) > 0x7e) => {
                let mut buf = [0u16; 2];
                for u in c.encode_utf16(&mut buf) {
                    write!(o, "\\u{:04x}", u).unwrap();
                }
            }
            c => o.push(c),
        }
    }
    o + "\""
}

fn jlist(items: impl Iterator<Item = String>) -> String {
    format!("[{}]", items.collect::<Vec<_>>().join(", "))
}

/// Every link but an in-page #anchor opens in a new tab: each
/// `<a\b[^>]*>` whose href\s*=\s*["']?([^"'\s>]*) is not #... or ?..., and
/// which has no \btarget\s*=.
fn new_tab(html: &str) -> String {
    let word = |c: char| c.is_alphanumeric() || c == '_';
    let mut out = String::new();
    let mut rest = html;
    while let Some(i) = rest.find("<a") {
        let after = rest[i + 2..].chars().next();
        let end = rest[i..].find('>');
        if after.is_some_and(word) || end.is_none() {
            out += &rest[..i + 2];
            rest = &rest[i + 2..];
            continue;
        }
        let end = i + end.unwrap();
        out += &rest[..i];
        out += &fix_tag(&rest[i..=end]);
        rest = &rest[end + 1..];
    }
    out + rest
}

fn fix_tag(tag: &str) -> String {
    let ws = |c: char| c.is_whitespace();
    let skip_ws = |s: &str| s.trim_start_matches(ws).len();
    let mut href = None;
    for (s, _) in tag.match_indices("href") {
        let r = &tag[s + 4..];
        let r = &r[r.len() - skip_ws(r)..];
        let Some(r) = r.strip_prefix('=') else { continue };
        let r = &r[r.len() - skip_ws(r)..];
        let r = r.strip_prefix(['"', '\'']).unwrap_or(r);
        let v: String = r.chars().take_while(|&c| c != '"' && c != '\'' && c != '>' && !ws(c)).collect();
        href = Some(v);
        break;
    }
    let has_target = tag.match_indices("target").any(|(s, _)| {
        let before = tag[..s].chars().next_back().is_none_or(|c| !(c.is_alphanumeric() || c == '_'));
        let r = &tag[s + 6..];
        before && r.trim_start_matches(ws).starts_with('=')
    });
    match href {
        Some(h) if !h.starts_with('#') && !h.starts_with('?') && !has_target => format!("{} target=\"_blank\" rel=\"noopener\">", &tag[..tag.len() - 1]),
        _ => tag.to_string(),
    }
}

fn page(kerf: f64) -> String {
    let ps = pieces(kerf);
    let piece = |n: &str| ps.iter().find(|p| p.name == n).unwrap();
    let track: String = ruler_track(300).iter().map(|d| d.to_string()).collect();
    let disc: String = de_bruijn(2, 5).iter().map(|d| d.to_string()).collect();
    let (sdb, s3, bits3) = stack_db_html(3);
    let grids = grid_rings();
    let (ladder, vcards, tess) = venn_ladder();
    let (even, odd) = grid_sides();
    let ring = grid_ring();
    let sub: Vec<(&str, String)> = vec![
        ("{{RULER}}", piece("window-ruler-300").preview("the window ruler")),
        ("{{READER}}", piece("window-reader").preview("the reader card")),
        ("{{GRAY_DISC}}", piece("gray-disc").preview("the Gray code disc")),
        ("{{DB_DISC}}", piece("de-bruijn-disc").preview("the de Bruijn disc")),
        ("{{TRACK}}", track),
        ("{{DISC}}", disc),
        ("{{COUNTING}}", counting_both(4)),
        ("{{STACK_DB}}", sdb),
        ("{{STACK_BITS}}", bits3),
        ("{{VENN_LADDER}}", ladder),
        ("{{VENN_CARDS}}", vcards),
        ("{{TESSERACT}}", tess),
        ("{{GRID_LIST}}", grid_list(&ring)),
        ("{{EVEN}}", even.join(" ")),
        ("{{ODD}}", odd.join(" ")),
        ("{{LOGIC_JSON}}", jlist(logic_rings().iter().map(|r| jlist(r.iter().map(|t| js(t, true)))))),
        ("{{GRID_JSON}}", jlist(ring.iter().map(|k| { let t = LETTERS[index_of(k)].1; format!("[{}, [{}, {}]]", js(k, true), js(t.0, true), js(t.1, true)) }))),
        ("{{STACK_JSON}}", jlist(s3.iter().map(|c| js(&c.to_string(), true)))),
        ("{{MEANING_JSON}}", format!("{{{}}}", MEANING.iter().map(|(b, m)| format!("{}: {}", js(b, true), js(m, true))).collect::<Vec<_>>().join(", "))),
        ("{{GRIDS_JSON}}", jlist(grids.iter().map(|r| jlist(r.iter().map(|t| js(t, true)))))),
        ("{{LETTERS_JSON}}", format!("{{{}}}", LETTERS.iter().map(|(k, (a, b))| format!("{}: [{}, {}]", js(k, false), js(a, false), js(b, false))).collect::<Vec<_>>().join(", "))),
        ("{{GRID_N}}", grids.len().to_string()),
        ("{{GRID_COUNT}}", GRID_COUNT.to_string()),
    ];
    let mut html = std::fs::read_to_string(here().join("page.html")).unwrap();
    for (a, b) in &sub {
        assert!(html.contains(a), "{a}");
        html = html.replace(a, b);
    }
    assert!(!html.contains("{{"));
    new_tab(&html)
}

fn die(msg: &str) -> ! {
    eprintln!("{msg}");
    std::process::exit(1)
}

fn main() {
    let (mut kerf, mut write) = (KERF, true);
    let dir = here();
    for a in std::env::args().skip(1) {
        if a.starts_with("--kerf=") {
            kerf = a.split('=').nth(1).unwrap().parse().unwrap_or_else(|_| die(&format!("could not convert string to float: '{}'", a.split('=').nth(1).unwrap())));
        } else if a == "--check" {
            write = false;
        } else if a == "--resample" {
            let rings = sample_grid_rings(48, 1);
            let body = rings.iter().map(|r| format!("[\n{}\n]", r.iter().map(|t| js(t, false)).collect::<Vec<_>>().join(",\n"))).collect::<Vec<_>>().join(",\n");
            std::fs::write(dir.join("grid-rings.json"), format!("[\n{body}\n]")).unwrap();
            println!("wrote grid-rings.json");
        } else if a == "--graph" {
            // the least-turning search's graph, for tools/count_grid_rings.rs
            let nb = grid_moves();
            println!("{} {}", LETTERS.len(), index_of("G"));
            for v in &nb {
                let mut v = v.clone();
                v.sort_by(|x, y| (letter(x.0), x.1).cmp(&(letter(y.0), y.1)));
                println!("{}", v.iter().filter(|e| e.1 <= 2).map(|e| format!("{}:{}", e.0, e.1)).collect::<Vec<_>>().join(" "));
            }
            return;
        } else {
            die(&format!("unknown {a}; takes --kerf=MM, --check, --resample, --graph"));
        }
    }
    if !check() {
        die("checks failed");
    }
    for p in pieces(kerf) {
        let (w, h) = p.size();
        let r3 = |v: f64| (v * 1000.0).round() / 1000.0;
        let fits = r3(w) <= BED_W && r3(h) <= BED_H;
        println!("{:<20} {w:7.2} x {h:6.2} mm  {} the P2S bed", p.name, if fits { "fits" } else { "DOES NOT FIT" });
        if !fits {
            std::process::exit(1);
        }
        if write {
            std::fs::create_dir_all(dir.join("cut-files")).unwrap();
            std::fs::write(Path::new(&dir).join("cut-files").join(format!("{}.svg", p.name)), p.cut_svg()).unwrap();
        }
    }
    if write {
        std::fs::write(dir.join("index.html"), page(KERF)).unwrap();
        println!("wrote cut-files/ (kerf {}) and index.html", py_float(kerf));
    }
}
