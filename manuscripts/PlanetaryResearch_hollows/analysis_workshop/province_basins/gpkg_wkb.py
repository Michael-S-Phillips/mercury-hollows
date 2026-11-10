import sqlite3, struct, numpy as np
def parse_gpb(b):
    flags = b[3]; env = (flags >> 1) & 7; n_env = {0:0,1:4,2:6,3:6,4:8}[env]
    off = 8 + 8*n_env; return parse_wkb(b, off)[0]
def parse_wkb(b, off):
    bo = '<' if b[off]==1 else '>'; t = struct.unpack(bo+'I', b[off+1:off+5])[0]; off += 5
    dim = 3 if (t//1000 in (1,2)) else (4 if t//1000==3 else 2); base = t % 1000
    if base == 2:
        n = struct.unpack(bo+'I', b[off:off+4])[0]; off += 4
        a = np.frombuffer(b[off:off+8*dim*n], dtype=bo+'f8').reshape(n, dim)[:, :2]; return [a], off + 8*dim*n
    if base == 5:
        n = struct.unpack(bo+'I', b[off:off+4])[0]; off += 4; out = []
        for _ in range(n):
            r, off = parse_wkb(b, off); out += r
        return out, off
    raise ValueError(t)
def load(path):
    c = sqlite3.connect(path); rows = c.execute("select fid, fault_type, sph_len_m, geom from tectonic_lines").fetchall()
    return [(fid, ft, L, parse_gpb(g)) for fid, ft, L, g in rows]
