import struct, pandas as pd
def read_dbf(path):
    with open(path, "rb") as f:
        hdr = f.read(32); n = struct.unpack("<I", hdr[4:8])[0]; hl, rl = struct.unpack("<HH", hdr[8:12])
        fields = []
        while True:
            fd = f.read(32)
            if fd[0] == 0x0D: break
            fields.append((fd[:11].split(b"\0")[0].decode(), chr(fd[11]), fd[16]))
        f.seek(hl); rows = []
        for _ in range(n):
            r = f.read(rl); off = 1; row = {}
            for name, t, L in fields:
                v = r[off:off + L].decode("latin1").strip(); off += L
                if t in "NF":
                    try: v = float(v)
                    except: v = None
                row[name] = v
            rows.append(row)
    return pd.DataFrame(rows)
