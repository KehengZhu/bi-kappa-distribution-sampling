"""Binary and JSONL readers for Experiment 6 raw output.

Every bulk file carries a header naming its schema version and record size.  The readers
refuse a file whose header does not match what they were written to understand, and
cross-check the declared record count against the file length, so a truncated or
mis-specified file fails loudly instead of being parsed into plausible nonsense.
"""

from __future__ import annotations

import json
import os
import struct

import numpy as np

HEADER_FMT = "<8sIIII7dIIq"
HEADER_SIZE = struct.calcsize(HEADER_FMT)  # 96

RECORD_KINDS = {1: "pilot", 2: "conditioning", 3: "loader", 4: "audit"}

PILOT_DTYPE = np.dtype([("log_r", "<f8"), ("log_w", "<f8"),
                        ("retries", "<u4"), ("flags", "<u4")])
COND_DTYPE = np.dtype([("log_r_ref", "<f8"), ("log_w_ref", "<f8"),
                       ("log_speed_ref", "<f8"), ("max_log_component_ref", "<f8"),
                       ("flags", "<u4"), ("cat_qf", "u1"), ("cat_split", "u1"),
                       ("cat_log", "u1"), ("reserved", "u1")])
LOADER_DTYPE = np.dtype([("v", "<f8", 3), ("log_r_ref", "<f8"),
                         ("status", "<u4"), ("attempts", "<u4")])
AUDIT_DTYPE = np.dtype([("x1", "<f8"), ("y", "<f8"), ("u", "<f8"), ("cos_theta", "<f8"),
                        ("phi", "<f8"), ("log_x2_working", "<f8"), ("kappa", "<f8"),
                        ("flags", "<u4"), ("precision_is_float", "<u4"),
                        ("cat_qf", "u1"), ("cat_split", "u1"), ("cat_log", "u1"),
                        ("reserved", "u1"), ("pad", "<u4")])

DTYPE_BY_KIND = {1: PILOT_DTYPE, 2: COND_DTYPE, 3: LOADER_DTYPE, 4: AUDIT_DTYPE}

CATEGORY_NAMES = ["finite", "denominator_zero", "quotient_first_loss", "split_form_loss",
                  "log_primitive_failure", "honest_overflow", "cap_reject", "cap_exhausted"]

FLAG_BITS = {
    "x2_zero": 1 << 0, "x2_subnormal": 1 << 1, "qf_nonfinite": 1 << 2,
    "split_nonfinite": 1 << 3, "log_nonfinite": 1 << 4, "honest_overflow_ref": 1 << 5,
    "cap_accept": 1 << 6, "audited": 1 << 7, "near_limit": 1 << 8,
    "u_endpoint_redraw": 1 << 9, "rotation_recoverable": 1 << 10, "cap_exhausted": 1 << 11,
}


class SchemaError(RuntimeError):
    pass


def read_records(path: str, expect_kind: int, expect_schema: int = 1):
    """Return ``(header_dict, structured_array)`` for one bulk file."""
    size = os.path.getsize(path)
    with open(path, "rb") as fh:
        raw = fh.read(HEADER_SIZE)
        if len(raw) < HEADER_SIZE:
            raise SchemaError(f"{path}: shorter than one header")
        (magic, schema, kind, rsize, _res, kappa, tperp, tpar, ub0, ub1, ub2, cap,
         seed, is_float, n_records) = struct.unpack(HEADER_FMT, raw)
        if magic[:7] != b"EXP6REC":
            raise SchemaError(f"{path}: not an exp6 record file")
        if kind != expect_kind:
            raise SchemaError(f"{path}: record kind {kind}, expected {expect_kind}")
        if schema != expect_schema:
            raise SchemaError(f"{path}: schema {schema}, expected {expect_schema}")
        dt = DTYPE_BY_KIND[kind]
        if rsize != dt.itemsize:
            raise SchemaError(f"{path}: record size {rsize}, reader expects {dt.itemsize}")
        expected_bytes = HEADER_SIZE + n_records * rsize
        if expected_bytes != size:
            raise SchemaError(
                f"{path}: header declares {n_records} records ({expected_bytes} bytes) "
                f"but the file is {size} bytes -- truncated or still being written")
        arr = np.fromfile(fh, dtype=dt, count=n_records)
    header = {"schema": schema, "kind": RECORD_KINDS[kind], "kappa": kappa,
              "theta_perp": tperp, "theta_par": tpar, "ub": (ub0, ub1, ub2), "cap": cap,
              "seed": seed, "precision": "float" if is_float else "double",
              "n_records": n_records, "path": path}
    return header, arr


def read_jsonl(path: str) -> list[dict]:
    rows = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def load_phase(raw_dir: str, phase: str, smoke: bool = False) -> list[dict]:
    """Every JSONL counter row a phase produced, across toolchains.

    Smoke output lives in its own directory and is never mixed in unless it is asked for
    explicitly: a smoke row carries 1000 attempts and would otherwise sit in a pooled rate
    beside a production row carrying a million.
    """
    sub = "smoke" if smoke else phase
    d = os.path.join(raw_dir, sub)
    if not os.path.isdir(d):
        return []
    rows = []
    for name in sorted(os.listdir(d)):
        if name.startswith(phase + "_") and name.endswith(".jsonl"):
            rows.extend(read_jsonl(os.path.join(d, name)))
    return rows


def field_basis(ub) -> np.ndarray:
    """The field-aligned basis exactly as bi_kappa_distribution builds it.

    Re-derived here rather than imported, so that a frame test is a genuine independent
    check of the rotation rather than a tautology.  Columns are (e1, e2, e3); the released
    header returns ``local[0] e1 + local[1] e2 + local[2] e3``.
    """
    ub = np.asarray(ub, dtype=float)
    e3 = ub / np.sqrt(ub @ ub)
    maxcomp = int(np.argmax(np.abs(e3)))
    e2 = np.ones(3)
    e2[maxcomp] = 1.0 - e3.sum() / e3[maxcomp]
    e2 = e2 / np.sqrt(e2 @ e2)
    e1 = np.cross(e2, e3)
    return np.column_stack([e1, e2, e3])


def recover_radius_direction(v: np.ndarray, kappa: float, theta_perp: float,
                             theta_par: float, ub) -> tuple[np.ndarray, np.ndarray]:
    """Recover ``(log R, unit direction in the field-aligned frame)`` from returned vectors.

    Done through a per-sample rescaling so that nothing overflows: the returned components
    reach 1e308 at the lowest kappa, and forming ``Q^T v`` directly would overflow on draws
    the loader successfully returned.
    """
    v = np.asarray(v, dtype=float)
    q = field_basis(ub)
    m = np.max(np.abs(v), axis=1)
    ok = np.isfinite(m) & (m > 0)
    log_r = np.full(v.shape[0], np.nan)
    n_hat = np.full(v.shape, np.nan)
    if not np.any(ok):
        return log_r, n_hat
    u = v[ok] / m[ok, None]              # order unity
    local = u @ q                        # = Q^T u, columns of q are the basis vectors
    scale = np.array([np.sqrt(kappa) * theta_perp, np.sqrt(kappa) * theta_perp,
                      np.sqrt(kappa) * theta_par])
    d = local / scale                    # order unity, direction times R/m
    s = np.hypot(np.hypot(d[:, 0], d[:, 1]), d[:, 2])
    good = s > 0
    idx = np.where(ok)[0][good]
    log_r[idx] = np.log(m[ok][good]) + np.log(s[good])
    n_hat[idx] = d[good] / s[good, None]
    return log_r, n_hat
