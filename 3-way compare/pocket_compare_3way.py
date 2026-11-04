# pocket_compare_3way.py -> PyMOL: cd C:/Users/NELSON/Desktop  then  run pocket_compare_3way.py
# Needs 3 objects loaded: 3T0H-empty, 5CF0-clean (ligand FJS), 3K60-clean (ligand ADP)
from pymol import cmd
import math, csv, os, re
from collections import OrderedDict

# ===================== USER SETTINGS =====================
NAMES   = {"t0h": "3T0H-empty", "h5cf0": "5CF0-clean", "k60": "3K60-clean"}  # short name : name as loaded
CHAIN   = "A"
LIG_H   = "resn FJS"   # ligand in 5CF0
LIG_K   = "resn ADP"   # ligand in 3K60
CUTOFF  = 5.0          # A, pocket = residues with ANY atom within this of the ligand
MOBILE_D = 1.5         # A, human side-chain shift above this = mobile in humans
SHIFT_D  = 2.5         # A, identical residue but 3K60 side chain this far from 3T0H = SHIFTED
FLAG_D   = 3.5         # A, CA-CA above this is flagged (loop moved)
OWN_D    = 2.5         # A, 3K60 further than this from BOTH human structures = own conformation
MARGIN_D = 1.0         # A, 3K60 must be this much closer to one human structure to be called like it
CONTACT_D = 4.5        # A, ligand atom this close to a residue = contact
EDGE_D    = 6.0        # A, up to this = edge (reachable by growing the ligand)
# Manual register fixes: (first 3T0H residue, last 3T0H residue, offset); 3T0H residue n pairs with 3K60 residue n-offset
# 132-138 (GQFGVGF) = 118-124 in 3K60, confirmed by identical sequence. CE had mis-paired Gly121 with Gly132.
REGISTER = [(132, 138, 14)]
OUTDIR  = os.getcwd()
# =========================================================

AA = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
      'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
KD  = dict(A=1.8,R=-4.5,N=-3.5,D=-3.5,C=2.5,Q=-3.5,E=-3.5,G=-0.4,H=-3.2,I=4.5,L=3.8,K=-3.9,
           M=1.9,F=2.8,P=-1.6,S=-0.8,T=-0.7,W=-0.9,Y=-1.3,V=4.2)
VOL = dict(A=88.6,R=173.4,N=114.1,D=111.1,C=108.5,Q=143.8,E=138.4,G=60.1,H=153.2,I=166.7,L=166.7,
           K=168.6,M=162.9,F=189.9,P=112.7,S=89.0,T=116.1,W=227.8,Y=193.6,V=140.0)
DON = dict(R=3,K=1,N=1,Q=1,H=1,S=1,T=1,Y=1,W=1)
ACC = dict(D=2,E=2,N=1,Q=1,H=1,S=1,T=1,Y=1)
CHG = dict(D=-1,E=-1,K=1,R=1)
AROM = set("FWYH")
GROUPS = ["AVILM","FWY","ST","NQ","DE","KRH","G","P","C"]
BACKBONE = {"N","CA","C","O","OXT"}

def grp(a): return next((g for g in GROUPS if a in g), a)
def pclass(a):
    if a == "G": return "Gly"
    if a == "P": return "Pro"
    if a in "DE": return "acidic(-)"
    if a in "KR": return "basic(+)"
    if a == "H": return "His(+/0)"
    if a in "AVILMFWC": return "hydrophobic"
    return "polar"
def ionis(a):
    return {"D":"acid","E":"acid","K":"base","R":"base","H":"His(pH)","Y":"weak(Tyr)","C":"weak(Cys)"}.get(a,"-")
def category(x, y):
    if x is None or y is None: return "NO EQUIVALENT"
    if x == y: return "IDENTICAL"
    return "SIMILAR" if grp(x) == grp(y) else "DIFFERENT"
def change(a, b):
    """what changes going from the HUMAN residue a to the 3K60 residue b"""
    if a is None or b is None: return "no equivalent residue"
    if a == b: return "-"
    s = []
    ca, cb = CHG.get(a,0), CHG.get(b,0)
    if ca != cb:
        if ca*cb < 0: s.append("charge reversal (%+d->%+d)" % (ca,cb))
        elif cb == 0: s.append("loses %s charge" % ("negative" if ca < 0 else "positive"))
        else:         s.append("gains %s charge" % ("negative" if cb < 0 else "positive"))
    d = DON.get(b,0) - DON.get(a,0)
    if d: s.append("%+d donor" % d)
    d = ACC.get(b,0) - ACC.get(a,0)
    if d: s.append("%+d acceptor" % d)
    if (a in AROM) != (b in AROM): s.append("gains aromatic ring" if b in AROM else "loses aromatic ring")
    dk = KD[b] - KD[a]
    if abs(dk) >= 3: s.append("more hydrophobic" if dk > 0 else "more polar")
    dv = VOL[b] - VOL[a]
    if abs(dv) >= 20: s.append("bulkier (%+.0f A3)" % dv if dv > 0 else "smaller (%+.0f A3)" % dv)
    return "; ".join(s) if s else "minor change"
def numkey(resi):
    m = re.match(r'-?\d+', resi); return int(m.group()) if m else 0
def dist(p, q): return math.dist(p, q) if p and q else None
def lab(r): return "%s%s" % (r["resn"], r["resi"]) if r else "---"
def f2(v): return "%.2f" % v if v is not None else "-"
def f1(v): return "%.1f" % v if v is not None else "-"
def pr(aa):
    if aa is None: return ("-",)*6
    return ("%.1f" % KD[aa], pclass(aa), "Y" if aa in AROM else "-", ionis(aa), str(DON.get(aa,0)), str(ACC.get(aa,0)))
def mind(rec, xyz):
    """closest distance between any heavy atom of a residue and any ligand atom"""
    if not rec or not xyz: return None
    return min(math.dist(a, b) for a in rec["all"] for b in xyz)

def load(obj):
    m = cmd.get_model("%s and chain %s and polymer.protein" % (obj, CHAIN))
    d = OrderedDict()
    for at in m.atom:
        k = (at.chain, at.resi)
        if k not in d:
            d[k] = dict(obj=obj, key=k, resn=at.resn, resi=at.resi, aa=AA.get(at.resn,"X"),
                        ca=None, idx=None, sc=[], all=[])
        r = d[k]
        if at.symbol != "H": r["all"].append(tuple(at.coord))
        if at.name == "CA" and r["ca"] is None:
            r["ca"] = tuple(at.coord); r["idx"] = at.index
        elif at.name not in BACKBONE and at.symbol != "H":
            r["sc"].append(at.coord)
    out = OrderedDict()
    for k, r in d.items():
        if r["ca"] is None: continue
        pts = r["sc"]
        r["cen"] = tuple(sum(p[i] for p in pts)/len(pts) for i in range(3)) if pts else r["ca"]
        out[k] = r
    return out

def pairs(aln, o_t, o_m, res):
    i2k = {o: {r["idx"]: k for k, r in res[o].items()} for o in (o_t, o_m)}
    fwd, rev = {}, {}
    for col in cmd.get_raw_alignment(aln):
        d = dict(col)
        if o_t in d and o_m in d:
            a, b = i2k[o_t].get(d[o_t]), i2k[o_m].get(d[o_m])
            if a and b: fwd[a] = b; rev[b] = a
    return fwd, rev

def main():
    names = cmd.get_names("objects")
    for new, old in NAMES.items():
        if new not in names and old in names: cmd.set_name(old, new)
    names = cmd.get_names("objects")
    miss = [n for n in NAMES if n not in names]
    if miss:
        print("ERROR: objects not found: %s. Loaded: %s" % (miss, names)); return
    L = []
    def out(s=""): print(s); L.append(s)
    cmd.delete("aln_k aln_h probe_fjs probe_adp pair_lines tmp_pk")
    cmd.label("all", "")

    sel = lambda o: "%s and chain %s and polymer.protein" % (o, CHAIN)
    # everything is superimposed onto 3T0H, so 3T0H numbering is the common key
    r1 = cmd.cealign(sel("t0h"), sel("k60"),   object="aln_k")
    r2 = cmd.cealign(sel("t0h"), sel("h5cf0"), object="aln_h")
    out("CE alignment 3K60 -> 3T0H : RMSD %.2f A over %d residues" % (r1["RMSD"], r1["alignment_length"]))
    out("CE alignment 5CF0 -> 3T0H : RMSD %.2f A over %d residues" % (r2["RMSD"], r2["alignment_length"]))

    res = {o: load(o) for o in NAMES}
    fk, rk = pairs("aln_k", "t0h", "k60", res)       # 3T0H->3K60, 3K60->3T0H
    fh, rh = pairs("aln_h", "t0h", "h5cf0", res)     # 3T0H->5CF0, 5CF0->3T0H

    # manual register fixes (override the CE pairing)
    for lo, hi, off in REGISTER:
        for n in range(lo, hi + 1):
            tk, kk = (CHAIN, str(n)), (CHAIN, str(n - off))
            if tk in res["t0h"] and kk in res["k60"]:
                oldk = fk.pop(tk, None)
                if oldk: rk.pop(oldk, None)
                oldt = rk.pop(kk, None)
                if oldt: fk.pop(oldt, None)
                fk[tk] = kk; rk[kk] = tk
        out("Register override applied: 3T0H %d-%d paired with 3K60 %d-%d (offset %d)" % (lo, hi, lo-off, hi-off, off))

    cmd.select("probe_fjs", "h5cf0 and chain %s and %s" % (CHAIN, LIG_H))
    cmd.select("probe_adp", "k60 and chain %s and %s" % (CHAIN, LIG_K))
    if cmd.count_atoms("probe_fjs") == 0 or cmd.count_atoms("probe_adp") == 0:
        print("ERROR: empty ligand selection (FJS atoms: %d, ADP atoms: %d)" %
              (cmd.count_atoms("probe_fjs"), cmd.count_atoms("probe_adp"))); return
    fjs_xyz = [tuple(a.coord) for a in cmd.get_model("probe_fjs").atom]
    adp_xyz = [tuple(a.coord) for a in cmd.get_model("probe_adp").atom]

    def pocket(o, probe):
        # select whole residues with ANY atom within CUTOFF, then take their CA atoms
        cmd.select("tmp_pk", "byres (%s within %s of %s)" % (sel(o), CUTOFF, probe))
        m = cmd.get_model("tmp_pk and name CA")
        cmd.delete("tmp_pk")
        return {(a.chain, a.resi) for a in m.atom} & set(res[o])
    pk = {"t_f": pocket("t0h", "probe_fjs"), "h_f": pocket("h5cf0", "probe_fjs"),
          "t_a": pocket("t0h", "probe_adp"), "k_a": pocket("k60", "probe_adp")}

    T = set(pk["t_f"]) | set(pk["t_a"])
    T |= {rh[k] for k in pk["h_f"] if k in rh}
    T |= {rk[k] for k in pk["k_a"] if k in rk}
    orphans = sorted([k for k in pk["k_a"] if k not in rk], key=lambda k: numkey(k[1]))

    pos = []
    for t in sorted(T, key=lambda k: numkey(k[1])):
        rt = res["t0h"][t]
        h = res["h5cf0"].get(fh.get(t)); k = res["k60"].get(fk.get(t))
        fl = []
        if t in pk["t_f"]: fl.append("FJS:3T0H")
        if h and fh[t] in pk["h_f"]: fl.append("FJS:5CF0")
        if t in pk["t_a"]: fl.append("ADP:3T0H")
        if k and fk[t] in pk["k_a"]: fl.append("ADP:3K60")
        in_fjs = ("FJS:3T0H" in fl) or ("FJS:5CF0" in fl)
        dF_h, dF_t, dA_k = mind(h, fjs_xyz), mind(rt, fjs_xyz), mind(k, adp_xyz)
        dF = min([v for v in (dF_h, dF_t) if v is not None], default=None)
        reach = "far" if dF is None or dF > EDGE_D else ("contact" if dF <= CONTACT_D else "edge")
        d_hs = dist(rt["cen"], h["cen"]) if h else None
        d_tk_ca, d_tk_sc = (dist(rt["ca"], k["ca"]), dist(rt["cen"], k["cen"])) if k else (None, None)
        d_hk_ca = dist(h["ca"], k["ca"]) if (h and k) else None
        mobile = d_hs is not None and d_hs > MOBILE_D
        cat = category(rt["aa"], k["aa"] if k else None)
        if cat == "NO EQUIVALENT": st = "NO-3K60-PARTNER"
        elif cat == "DIFFERENT": st = "HANDLE-DIFFERENT"
        elif cat == "SIMILAR": st = "HANDLE-SIMILAR"
        elif d_tk_sc is not None and d_tk_sc > SHIFT_D: st = "SHIFTED"
        else: st = "CONSERVED"
        near = "-"
        if d_tk_ca is not None and d_hk_ca is not None:
            nearer = "3T0H" if d_tk_ca < d_hk_ca else "5CF0"
            if min(d_tk_ca, d_hk_ca) > OWN_D: near = "own conformation (nearer %s)" % nearer
            elif abs(d_tk_ca - d_hk_ca) >= MARGIN_D: near = "%s-like" % nearer
            else: near = "equal"
        pos.append(dict(t=rt, h=h, k=k, fl=" ".join(fl), in_fjs=in_fjs, d_hs=d_hs, d_tk_ca=d_tk_ca,
                        d_tk_sc=d_tk_sc, d_hk_ca=d_hk_ca, mobile=mobile, cat=cat, st=st, near=near,
                        x=rt["aa"], y=k["aa"] if k else None, dF_h=dF_h, dF_t=dF_t, dA_k=dA_k, reach=reach))
    for i, p in enumerate(pos, 1): p["n"] = i

    out("\nPocket = residues with ANY atom within %.1f A of FJS (in 5CF0) or ADP (in 3K60); all structures superimposed on 3T0H." % CUTOFF)
    out("Positions = union of all four pocket lists, keyed by 3T0H number. Change is read HUMAN -> 3K60.")
    W = "=" * 130
    out("\n" + W + "\nTABLE 1. POCKET RESIDUES\n" + W)
    for lbl, key, o in (("5CF0 (FJS)", "h_f", "h5cf0"), ("3T0H (FJS footprint)", "t_f", "t0h"),
                        ("3T0H (ADP footprint)", "t_a", "t0h"), ("3K60 (ADP)", "k_a", "k60")):
        ks = sorted(pk[key], key=lambda k: numkey(k[1]))
        out("%-22s (%2d): %s" % (lbl, len(ks), " ".join(lab(res[o][k]) for k in ks)))

    out("\n" + W + "\nTABLE 2. SAME POSITION IN 3D: 3T0H / 5CF0 / 3K60\n" + W)
    out("%-4s %-8s %-8s %-8s %-38s %-7s %-7s %-7s %-7s %s" %
        ("Pos", "3T0H", "5CF0", "3K60", "In pocket of", "T-H sc", "T-K CA", "T-K sc", "H-K CA", "3K60 loop looks like"))
    for p in pos:
        flag = "  [CHECK: loop moved]" if (p["d_tk_ca"] or 0) > FLAG_D else ""
        out("%-4s %-8s %-8s %-8s %-38s %-7s %-7s %-7s %-7s %s%s" %
            ("P%d" % p["n"], lab(p["t"]), lab(p["h"]), lab(p["k"]), p["fl"], f2(p["d_hs"]),
             f2(p["d_tk_ca"]), f2(p["d_tk_sc"]), f2(p["d_hk_ca"]), p["near"], flag))
    out("T-H sc = side-chain shift between the two human structures; T-K = 3T0H vs 3K60; H-K = 5CF0 vs 3K60 (A).")
    if orphans:
        out("3K60 ADP-pocket residues with no CE partner in 3T0H: " + " ".join(lab(res["k60"][k]) for k in orphans))

    out("\n" + W + "\nTABLE 3. SELECTIVITY MAP (human -> 3K60)\n" + W)
    out("%-4s %-8s %-8s %-18s %-7s %-8s %-8s %-8s %-8s %s" %
        ("Pos", "Human", "3K60", "Status", "Human", "FJS-d", "FJS-d", "ADP-d", "Reach", "Change going human -> 3K60"))
    out("%-4s %-8s %-8s %-18s %-7s %-8s %-8s %-8s %-8s %s" %
        ("", "", "", "", "", "(5CF0)", "(3T0H)", "(3K60)", "(FJS)", ""))
    for p in pos:
        out("%-4s %-8s %-8s %-18s %-7s %-8s %-8s %-8s %-8s %s" % ("P%d" % p["n"], lab(p["t"]), lab(p["k"]), p["st"],
                                              "MOBILE" if p["mobile"] else "rigid",
                                              f1(p["dF_h"]), f1(p["dF_t"]), f1(p["dA_k"]), p["reach"],
                                              change(p["x"], p["y"])))
    out("FJS-d / ADP-d = closest heavy-atom distance (A) from the ligand to the residue. Reach: contact <= %.1f A, edge <= %.1f A (reachable by growing the ligand), far = beyond." % (CONTACT_D, EDGE_D))
    out("HANDLE-DIFFERENT / HANDLE-SIMILAR = sequence differs in 3K60; SHIFTED = same residue but side chain >%.1f A away in 3K60;" % SHIFT_D)
    out("CONSERVED = no handle; NO-3K60-PARTNER = no CE partner, check by eye. 'rigid' human positions are the safest to exploit.")
    for title, keys in (("Handles, DIFFERENT residue", ("HANDLE-DIFFERENT",)), ("Handles, SIMILAR residue", ("HANDLE-SIMILAR",)),
                        ("SHIFTED side chains (shape handles)", ("SHIFTED",)), ("No 3K60 partner (check)", ("NO-3K60-PARTNER",))):
        sel_p = [p for p in pos if p["st"] in keys]
        out("\n%s:" % title)
        for p in sel_p:
            out("  %-8s -> %-8s %-6s reach=%-8s FJS-d=%s  %s" % (lab(p["t"]), lab(p["k"]), "MOBILE" if p["mobile"] else "rigid",
                                                  p["reach"], f1(min([v for v in (p["dF_h"], p["dF_t"]) if v is not None], default=None)),
                                                  change(p["x"], p["y"])))
        if not sel_p: out("  none")
    cons = [lab(p["t"]) for p in pos if p["st"] == "CONSERVED"]
    out("\nConserved in all three (potency anchors, no selectivity):\n  " + " ".join(cons))

    out("\n" + W + "\nTABLE 4. PROPERTIES AT EACH POSITION (side chains only)\n" + W)
    h = "%-4s | %-8s %5s %-12s %-3s %-10s %-2s %-2s | %-8s %5s %-12s %-3s %-10s %-2s %-2s"
    out(h % ("Pos", "Human", "KD", "Class", "Ar", "Ion", "D", "A", "3K60", "KD", "Class", "Ar", "Ion", "D", "A"))
    for p in pos:
        out(h % (("P%d" % p["n"], lab(p["t"])) + pr(p["x"]) + (lab(p["k"]),) + pr(p["y"])))
    out("Backbone N-H donors and C=O acceptors are not counted. KD = Kyte-Doolittle. Ar = aromatic. Ion = ionisable. D/A = donors/acceptors.")

    txt = os.path.join(OUTDIR, "threeway_report.txt")
    with open(txt, "w", encoding="utf-8") as fh: fh.write("\n".join(L))
    csvf = os.path.join(OUTDIR, "threeway_positions.csv")
    with open(csvf, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["pos", "3T0H", "5CF0", "3K60", "in_pocket_of", "FJS_d_5CF0", "FJS_d_3T0H", "ADP_d_3K60", "reach",
                    "human_sc_shift", "T_K_CA", "T_K_sc", "H_K_CA", "k60_looks_like", "category", "status",
                    "human_mobility", "change_human_to_3K60",
                    "human_KD", "human_class", "human_arom", "human_ion", "human_don", "human_acc",
                    "k60_KD", "k60_class", "k60_arom", "k60_ion", "k60_don", "k60_acc"])
        for p in pos:
            w.writerow(["P%d" % p["n"], lab(p["t"]), lab(p["h"]), lab(p["k"]), p["fl"], f1(p["dF_h"]), f1(p["dF_t"]),
                        f1(p["dA_k"]), p["reach"], f2(p["d_hs"]), f2(p["d_tk_ca"]), f2(p["d_tk_sc"]), f2(p["d_hk_ca"]),
                        p["near"], p["cat"], p["st"], "MOBILE" if p["mobile"] else "rigid", change(p["x"], p["y"])]
                       + list(pr(p["x"])) + list(pr(p["y"])))
    print("\nWritten:\n  %s\n  %s" % (txt, csvf))

    # display
    col = {"CONSERVED": "green", "HANDLE-DIFFERENT": "red", "HANDLE-SIMILAR": "yellow",
           "SHIFTED": "orange", "NO-3K60-PARTNER": "gray50"}
    cmd.hide("everything", "all")
    cmd.show("cartoon", " or ".join(NAMES)); cmd.set("cartoon_transparency", 0.75)
    cmd.color("wheat", "t0h"); cmd.color("lightblue", "h5cf0"); cmd.color("palegreen", "k60")
    cmd.show("sticks", "probe_fjs or probe_adp")
    cmd.color("magenta", "probe_fjs and elem C"); cmd.color("hotpink", "probe_adp and elem C")
    cmd.set("label_size", 14)
    for p in pos:
        for rec in (p["t"], p["h"], p["k"]):
            if not rec: continue
            s = "byres (%s and index %d)" % (rec["obj"], rec["idx"])
            cmd.show("sticks", s); cmd.color(col[p["st"]], s + " and elem C")
        cmd.label("t0h and index %d" % p["t"]["idx"], repr("P%d %s|%s" % (p["n"], lab(p["t"]), lab(p["k"]))))
    cmd.zoom("probe_fjs or probe_adp", 12); cmd.deselect()
    print("Colours: green=conserved, red=handle (different residue), yellow=handle (similar), orange=shifted side chain, grey=no 3K60 partner")
    print("Labels: P<n> 3T0H residue | 3K60 residue")

main()