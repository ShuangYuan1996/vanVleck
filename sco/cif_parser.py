"""Lightweight CIF 1.1 parser (no third-party dependencies).

Supports several data_ blocks, loop_, single/double-quoted strings, semicolon text fields and # comments.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

_NUM_RE = re.compile(r"^([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)(?:\(\d+\))?$")


def cif_float(value) -> float | None:
    """'0.60011(8)' -> 0.60011; '?' or '.' -> None."""
    if value is None:
        return None
    s = str(value).strip()
    if s in ("", "?", "."):
        return None
    m = _NUM_RE.match(s)
    if m:
        return float(m.group(1))
    try:
        return float(s)
    except ValueError:
        return None


@dataclass
class CifBlock:
    name: str
    items: dict = field(default_factory=dict)   # tag (lower case) -> str
    loops: list = field(default_factory=list)   # [dict tag (lower case) -> list[str]]

    def get(self, *tags, default=None):
        """Try equivalent tags in order (old and new CIF dictionaries); return the first value found."""
        for t in tags:
            v = self.items.get(t.lower())
            if v is not None and v not in ("?", "."):
                return v
        for t in tags:   # the item may also appear as a one-row loop
            loop = self.find_loop(t)
            if loop and loop[t.lower()]:
                v = loop[t.lower()][0]
                if v not in ("?", "."):
                    return v
        return default

    def find_loop(self, *tags):
        for t in tags:
            t = t.lower()
            for loop in self.loops:
                if t in loop:
                    return loop
        return None

    def column(self, loop, *tags):
        for t in tags:
            if t.lower() in loop:
                return loop[t.lower()]
        return None


def _tokenize(text: str):
    """Yield (kind, value) with kind in DATA / LOOP / TAG / VALUE."""
    lines = text.splitlines()
    i, n_lines = 0, len(lines)
    while i < n_lines:
        line = lines[i]
        if line.startswith(";"):
            buf = [line[1:]]
            i += 1
            while i < n_lines and not lines[i].startswith(";"):
                buf.append(lines[i])
                i += 1
            i += 1
            yield ("VALUE", "\n".join(buf).strip())
            continue
        pos, n = 0, len(line)
        while pos < n:
            c = line[pos]
            if c.isspace():
                pos += 1
                continue
            if c == "#":
                break
            if c in "'\"":
                j = pos + 1
                while True:
                    k = line.find(c, j)
                    if k == -1:
                        val, pos = line[pos + 1:], n
                        break
                    if k + 1 >= n or line[k + 1].isspace():
                        val, pos = line[pos + 1:k], k + 1
                        break
                    j = k + 1
                yield ("VALUE", val)
                continue
            j = pos
            while j < n and not line[j].isspace():
                j += 1
            tok = line[pos:j]
            pos = j
            low = tok.lower()
            if low.startswith("data_"):
                yield ("DATA", tok[5:])
            elif low == "loop_":
                yield ("LOOP", None)
            elif low.startswith("save_") or low in ("global_", "stop_"):
                continue
            elif tok.startswith("_"):
                yield ("TAG", low)
            else:
                yield ("VALUE", tok)
        i += 1


def parse_cif(text: str) -> list[CifBlock]:
    tokens = list(_tokenize(text))
    blocks: list[CifBlock] = []
    cur: CifBlock | None = None
    idx, n = 0, len(tokens)
    while idx < n:
        kind, val = tokens[idx]
        if kind == "DATA":
            cur = CifBlock(name=val)
            blocks.append(cur)
            idx += 1
        elif cur is None:
            idx += 1
        elif kind == "LOOP":
            idx += 1
            tags = []
            while idx < n and tokens[idx][0] == "TAG":
                tags.append(tokens[idx][1])
                idx += 1
            values = []
            while idx < n and tokens[idx][0] == "VALUE":
                values.append(tokens[idx][1])
                idx += 1
            if tags:
                nrow = len(values) // len(tags)
                loop = {t: [values[r * len(tags) + c] for r in range(nrow)]
                        for c, t in enumerate(tags)}
                cur.loops.append(loop)
        elif kind == "TAG":
            if idx + 1 < n and tokens[idx + 1][0] == "VALUE":
                cur.items[val] = tokens[idx + 1][1]
                idx += 2
            else:
                cur.items[val] = None
                idx += 1
        else:
            idx += 1
    return blocks
