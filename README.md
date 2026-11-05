# Gray and de Bruijn

**[← All projects](https://gernreich.github.io/)**

**[Read it, with a ruler and a disc you can decode](https://gernreich.github.io/gray-and-de-bruijn/)**

Two classic ways to walk through a small alphabet, tried on three alphabets
from earlier projects:

- **A Gray code** lists every symbol once, each differing from the one before
  by a single change.
- **A de Bruijn sequence** is one cyclic string in which every word of *n*
  symbols turns up once, so a window *n* wide tells you where you are.

The alphabets are the A B C layers of
[atomic layer memory](https://github.com/Gernreich/atomic-layer-memory), the
sixteen truth tables of Shea Zellweger's Logic Alphabet (as in
[Logic-Research](https://github.com/Gernreich/Logic-Research)) and the 28
letters of the [grid script](https://github.com/Gernreich/gernreich-writing).

What comes out:

- **A window ruler.** 300 mm, with a track of 1 mm cells along one edge taken
  from the de Bruijn sequence for windows of nine. Any nine cells seen through
  the reader card name the millimetre, at all 301 positions.
- **Two encoder discs.** 32 angles each: one with five Gray-code rings, one
  with a single de Bruijn ring read five cells at a time.
- **Layer memory.** Counting 0 to 15 moves 26 layers in binary and 15 in Gray
  code. A stack in which every possible window of *n* layers turns up once is
  3 × 2<sup>*n*−1</sup> layers long. Read as turns up or down (Hägg notation),
  a stack holds one bit per layer, twice the density of the spacer code.
- **A Logic Alphabet ring.** There are 1344 Gray cycles through the sixteen
  truth tables (Hamiltonian cycles on the 4-cube). In 48 of them every glyph
  sits across the ring from its negation. The ring here is one of those 48.
- **A grid-script ring.** No Gray cycle through the 28 letters can turn one
  step by only an eighth each time. Such a move always switches between
  even-bend and odd-bend letters, and there are 12 of one and 16 of the other.
  The least turning possible is 32 eighths, and the ring reaches it.

Every count above is checked by `walks.py`.

| Path | What |
|---|---|
| `walks.py` | the walks, the checks, the cut files and the page: `python3 walks.py` (`--kerf=MM`, default 0.13; `--check` writes nothing) |
| `page.html` | the page template; `walks.py` writes `index.html` from it |
| `index.html` | the page, with two 3D views (three.js): the window stack as close-packed atoms and the Logic ring on a turning 4-cube |
| `cut-files/` | `window-ruler-300`, `window-reader`, `gray-disc`, `de-bruijn-disc`, `logic-ring`, `grid-ring`: SVG for the xTool P2S. Green engrave, blue score, orange inner cut, black outline |

Gray codes are named for Frank Gray (patented 1953), though Émile Baudot used
one in the 1870s. De Bruijn sequences are named for Nicolaas de Bruijn (1946),
and Camille Flye Sainte-Marie found them first (1894). The layer memory and the
grid script are the author's, and the Logic Alphabet is Zellweger's.

Released under [CC0 1.0](LICENSE).

## Credits

Made by Michael Gernreich with [Claude](https://www.anthropic.com/claude) (Anthropic), working in Claude Code. Claude proposed this exploration at Michael's request, worked out the results (the 48 antipodal Logic rings, the parity bound of 32 eighths on the grid script, the stack counts), designed the pieces, and wrote `walks.py`, the cut files, the 3D views and the documentation. Commits Claude made carry a `Co-Authored-By: Claude` line.
