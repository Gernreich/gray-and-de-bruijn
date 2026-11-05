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
- **Layer memory.** The B spacers are a buffer and are ignored throughout.
  Counting 0 to 15 moves 26 bit layers in binary and 15 in Gray code. Written in
  the spacer code, the de Bruijn sequence 00010111 makes a 16-layer stack in
  which every 3-bit word turns up once among the bit layers.
- **A Logic Alphabet ring.** There are 1344 Gray cycles through the sixteen
  truth tables (Hamiltonian cycles on the 4-cube). In 48 of them every glyph
  sits across the ring from its negation. The page shows all 48, with a table
  that spells out each step. Ring 1 is the one on the cut file.
- **A grid-script ring.** No Gray cycle through the 28 letters can turn one
  step by only an eighth each time. Such a move always switches between
  even-bend and odd-bend letters, and there are 12 of one and 16 of the other.
  The least turning possible is 32 eighths, and the ring reaches it.

Every count above is checked by `walks.py`.

| Path | What |
|---|---|
| `walks.py` | the walks, the checks, the cut files and the page: `python3 walks.py` (`--kerf=MM`, default 0.13; `--check` writes nothing) |
| `page.html` | the page template; `walks.py` writes `index.html` from it |
| `index.html` | the page, with three 3D views (three.js): the spacer-code stack as close-packed atoms, the cut rings as 3 mm ply, and the Logic ring on a turning 4-cube |
| `cut-files/` | `window-ruler-300`, `window-reader`, `gray-disc`, `de-bruijn-disc`, `logic-ring`, `grid-ring`: SVG for the xTool P2S. Green engrave, blue score, orange inner cut, black outline |

Gray codes are named for Frank Gray (patented 1953), though Émile Baudot used
one in the 1870s. De Bruijn sequences are named for Nicolaas de Bruijn (1946),
and Camille Flye Sainte-Marie found them first (1894). The layer memory and the
grid script are the author's, and the Logic Alphabet is Zellweger's.

Released under [CC0 1.0](LICENSE).

## Credits

Made by Michael Gernreich with [Claude](https://www.anthropic.com/claude) (Anthropic), working in Claude Code. Claude proposed this exploration at Michael's request, worked out the results (the 48 antipodal Logic rings, the parity bound of 32 eighths on the grid script, the stack counts), designed the pieces, and wrote `walks.py`, the cut files, the 3D views and the documentation. Commits Claude made carry a `Co-Authored-By: Claude` line.
