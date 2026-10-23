"""테스트용 최소 ABIF(.ab1) 작성. Biopython SeqIO 가 PBAS/PCON 으로 읽는다."""
from __future__ import annotations

import struct
from array import array
from pathlib import Path


def write_ab1(path, sequence: str, qualities: list[int], sample: str = "sample") -> Path:
    seq = sequence.encode("ascii")
    qual = bytes(max(0, min(93, int(q))) for q in qualities)
    if len(seq) != len(qual):
        raise ValueError("sequence and qualities differ in length")
    step = 12
    nscan = step * len(seq) + 30
    data = [array("h", [0] * nscan) for _ in range(4)]
    ch = {ord("G"): 0, ord("A"): 1, ord("T"): 2, ord("C"): 3}
    ploc = []
    for i, base in enumerate(seq):
        pos = 15 + i * step
        ploc.append(pos)
        bi = ch.get(base)
        if bi is None:
            continue
        for d in range(-2, 3):
            p = pos + d
            if 0 <= p < nscan:
                data[bi][p] = 1500 - abs(d) * 400

    def shorts(vals) -> bytes:
        return b"".join(struct.pack(">h", int(x)) for x in vals)

    smpl = sample.encode("latin1")
    tags = [
        (b"PBAS", 1, 2, 1, len(seq), seq),
        (b"PBAS", 2, 2, 1, len(seq), seq),
        (b"PCON", 1, 2, 1, len(qual), qual),
        (b"PCON", 2, 2, 1, len(qual), qual),
        (b"FWO_", 1, 2, 1, 4, b"GATC"),
        (b"SMPL", 1, 2, 1, len(smpl), smpl),
        (b"PLOC", 1, 4, 2, len(ploc), shorts(ploc)),
        (b"PLOC", 2, 4, 2, len(ploc), shorts(ploc)),
    ]
    for i, arr in enumerate(data):
        payload = shorts(arr)
        tags.append((b"DATA", 9 + i, 4, 2, nscan, payload))

    dir_fmt = ">4sI2H4I"
    dir_size = struct.calcsize(dir_fmt)
    n = len(tags)
    header_len = 4 + struct.calcsize(">H4sI2H3I")
    cursor = header_len + n * dir_size
    entries = []
    blobs = []
    for name, number, code, esize, count, payload in tags:
        if len(payload) <= 4:
            inline = payload + b"\x00" * (4 - len(payload))
            off = struct.unpack(">I", inline)[0]
            entries.append(struct.pack(dir_fmt, name, number, code, esize, count, len(payload), off, 0))
        else:
            entries.append(struct.pack(dir_fmt, name, number, code, esize, count, len(payload), cursor, 0))
            blobs.append(payload)
            cursor += len(payload)
    header = b"ABIF" + struct.pack(
        ">H4sI2H3I", 101, b"tdir", 1, 1023, dir_size, n, n * dir_size, header_len
    )
    path = Path(path)
    path.write_bytes(header + b"".join(entries) + b"".join(blobs))
    return path
