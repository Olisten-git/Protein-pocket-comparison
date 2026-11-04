# pocket_compare_v2.py   ->   PyMOL:  cd C:/Users/NELSON/Desktop   then   run pocket_compare_v2.py
from pymol import cmd
import math, csv, os, re
from collections import OrderedDict

# ===================== USER SETTINGS =====================
REF_OLD, MOB_OLD = "3K60-clean", "3T0H-empty"   # names as you loaded them
REF, MOB = "k60", "t0h"                          # short names used internally
LBL_REF, LBL_MOB = "3K60", "3T0H"                # names used in the tables
CHAIN        = "A"
PROBE        = "resn ADP"    # ligand in REF that defines the pocket
INCLUDE_SO4  = False         # True: sulfates also define the pocket
CUTOFF       = 5.0           # A, pocket = residues with any atom within this of the probe
FALLBACK_D   = 4.0           # A, max CA-CA for 3D-nearest pairing when CE has no partner
FLAG_D       = 3.5           # A, CA-CA above this is flagged CHECK
DRAW_LINES   = True          # dashed lines between paired CA atoms
OUTDIR       = os.getcwd()
# =========================================================

AA = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
      'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
KD  = dict(A=1.8,R=-4.5,N=-3.5,D=-3.5,C=2.5,Q=-3.5,E=-3.5,G=-0.4,H=-3.2,I=4.5,L=3.8,K=-3.9,
           M=1.9,F=2.8,P=-1.6,S=-0.8,T=-0.7,W=-0.9,Y=-1.3,V=4.2)
VOL = dict(A=88.6,R=173.4,N=114.1,D=111.1,C=108.5,Q=143.8,E=138.4,G=60.1,H=153.2,I=166.7,L=166.7,
           K=168.6,M=162.9,F=189.9,P=112.7,S=89.0,T=116.1,W=227.8,Y=193.6,V=140.0)   # A^3
DON = dict(R=3,K=1,N=1,Q=1,H=1,S=1,T=1,Y=1,W=1)        # side-chain H-bond donors
ACC = dict(D=2,E=2,N=1,Q=1,H=1,S=1,T=1,Y=1)            # side-chain H-bond acceptors
CHG = dict(D=-1,E=-1,K=1,R=1)                          # charge at pH 7 (His treated as neutral)
AROM = set("FWYH")
GROUPS = ["AVILM","FWY","ST","NQ","DE","KRH","G","P","C"]   # "similar" = same group
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
    """what changes going from REF residue a to MOB residue b (design-oriented)"""
    if a is None or b is None: return "no equivalent residue at this position"
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

def load(obj):
    m = cmd.get_model("%s and chain %s and polymer.protein" % (obj, CHAIN))
    d = OrderedDict()
    for at in m.atom:
        k = (at.chain, at.resi)
        if k not in d:
            d[k] = dict(obj=obj, key=k, resn=at.resn, resi=at.resi, aa=AA.get(at.resn,"X"),
                        ca=None, idx=None, sc=[])
        r = d[k]
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

def nearest(rec, res, used):
    best, bd = None, 1e9
    for k, r in res.items():
        if k in used: continue
        dd = dist(rec["ca"], r["ca"])
        if dd < bd: best, bd = k, dd
    return best if bd <= FALLBACK_D else None

def totals(aas):
    aas = [x for x in aas if x]; n = len(aas)
    if n == 0: return OrderedDict([("residues", 0)])
    return OrderedDict([
        ("residues", n),
        ("hydrophobic (A,V,I,L,M,F,W,C)", sum(x in "AVILMFWC" for x in aas)),
        ("polar uncharged (S,T,N,Q,Y)",   sum(x in "STNQY" for x in aas)),
        ("acidic (D,E)",                   sum(x in "DE" for x in aas)),
        ("basic (K,R)",                    sum(x in "KR" for x in aas)),
        ("His",                            aas.count("H")),
        ("Gly",                            aas.count("G")),
        ("Pro",                            aas.count("P")),
        ("aromatic (F,W,Y,H)",             sum(x in AROM for x in aas)),
        ("net charge pH7 (His=0)",         sum(CHG.get(x,0) for x in aas)),
        ("mean Kyte-Doolittle",            round(sum(KD.get(x,0) for x in aas)/n, 2)),
        ("side-chain H-bond donors",       sum(DON.get(x,0) for x in aas)),
        ("side-chain H-bond acceptors",    sum(ACC.get(x,0) for x in aas))])

def main():
    names = cmd.get_names("objects")
    for old, new in ((REF_OLD, REF), (MOB_OLD, MOB)):
        if new not in names and old in names: cmd.set_name(old, new)
    names = cmd.get_names("objects")
    if REF not in names or MOB not in names:
        print("ERROR: could not find objects %s / %s. Loaded: %s" % (REF, MOB, names)); return

    L = []
    def out(s=""): print(s); L.append(s)

    # ---- clean up anything from earlier runs
    cmd.delete("aln pair_lines probe pocket_k60 pocket_t0h")
    cmd.label("all", ""); 

    # ---- 1. superimpose (CE alignment, which also gives the 1:1 residue pairing)
    tgt = "%s and chain %s and polymer.protein" % (REF, CHAIN)
    mob = "%s and chain %s and polymer.protein" % (MOB, CHAIN)
    r = cmd.cealign(tgt, mob, object="aln")
    out("CE alignment RMSD: %.2f A over %d residues" % (r["RMSD"], r["alignment_length"]))

    res = {REF: load(REF), MOB: load(MOB)}
    idx2key = {o: {r_["idx"]: k for k, r_ in res[o].items()} for o in res}
    ref2mob, mob2ref = {}, {}
    for col in cmd.get_raw_alignment("aln"):
        dd = {o: i for (o, i) in col}
        if REF in dd and MOB in dd:
            a, b = idx2key[REF].get(dd[REF]), idx2key[MOB].get(dd[MOB])
            if a and b: ref2mob[a] = b; mob2ref[b] = a
    out("Residues paired by the CE alignment: %d" % len(ref2mob))

    # ---- 2. pockets
    cmd.select("probe", "%s and chain %s and (%s%s)" % (REF, CHAIN, PROBE, " or resn SO4" if INCLUDE_SO4 else ""))
    if cmd.count_atoms("probe") == 0:
        print("ERROR: probe selection is empty"); return
    cmd.select("pocket_k60", "byres (%s within %s of probe)" % (tgt, CUTOFF))
    cmd.select("pocket_t0h", "byres (%s within %s of probe)" % (mob, CUTOFF))
    pk = {}
    for o, s in ((REF, "pocket_k60"), (MOB, "pocket_t0h")):
        keys = [(a.chain, a.resi) for a in cmd.get_model(s + " and name CA").atom]
        pk[o] = sorted([k for k in keys if k in res[o]], key=lambda k: numkey(k[1]))

    # ---- 3. build position table (same place in 3D)
    pos, used_a, used_b = [], set(), set()
    def add(a, b, how):
        ra, rb = res[REF].get(a) if a else None, res[MOB].get(b) if b else None
        if a: used_a.add(a)
        if b: used_b.add(b)
        inA, inB = bool(a and a in pk[REF]), bool(b and b in pk[MOB])
        where = "both pockets" if inA and inB else ("%s pocket only" % LBL_REF if inA else
                ("%s pocket only" % LBL_MOB if inB else "-"))
        dca = dist(ra["ca"], rb["ca"]) if ra and rb else None
        dsc = dist(ra["cen"], rb["cen"]) if ra and rb else None
        if dca is not None and dca > FLAG_D: how += "  [CHECK: CA-CA > %.1f A]" % FLAG_D
        pos.append(dict(a=ra, b=rb, where=where, how=how, dca=dca, dsc=dsc,
                        x=ra["aa"] if ra else None, y=rb["aa"] if rb else None))
    for a in pk[REF]:
        b = ref2mob.get(a); how = "CE alignment"
        if b is None or b in used_b:
            b = nearest(res[REF][a], res[MOB], used_b); how = "3D nearest CA" if b else "no partner"
        add(a, b, how)
    for b in pk[MOB]:
        if b in used_b: continue
        a = mob2ref.get(b); how = "CE alignment"
        if a is None or a in used_a:
            a = nearest(res[MOB][b], res[REF], used_a); how = "3D nearest CA" if a else "no partner"
        add(a, b, how)
    pos.sort(key=lambda p: numkey(p["a"]["resi"]) if p["a"] else 10000 + numkey(p["b"]["resi"]))
    for i, p in enumerate(pos, 1): p["n"] = i; p["cat"] = category(p["x"], p["y"])

    # ---- 4. report
    out("\nPocket = residues with any atom within %.1f A of %s (probe taken from %s, MOB superimposed on it)" %
        (CUTOFF, "ADP" + (" + SO4" if INCLUDE_SO4 else ""), LBL_REF))
    out("\n" + "="*100 + "\nTABLE 1. RESIDUES IN EACH POCKET\n" + "="*100)
    for lbl, o in ((LBL_REF, REF), (LBL_MOB, MOB)):
        out("%s (%d residues): %s" % (lbl, len(pk[o]), " ".join(lab(res[o][k]) for k in pk[o])))

    out("\n" + "="*100 + "\nTABLE 2. RESIDUES AT THE SAME POSITION IN 3D (aligned pockets)\n" + "="*100)
    out("%-4s %-9s %-9s %-22s %-9s %-14s %s" % ("Pos", LBL_REF, LBL_MOB, "In pocket of", "CA-CA(A)", "SideChain(A)", "Paired by"))
    for p in pos:
        f = lambda v: "%.2f" % v if v is not None else "-"
        out("%-4s %-9s %-9s %-22s %-9s %-14s %s" % ("P%d" % p["n"], lab(p["a"]), lab(p["b"]), p["where"],
                                                    f(p["dca"]), f(p["dsc"]), p["how"]))
    out("(SideChain = distance between side-chain centroids; a large value means the same sequence position but the side chain points elsewhere.)")

    out("\n" + "="*100 + "\nTABLE 3. IDENTICAL / SIMILAR / DIFFERENT AT EACH POSITION\n" + "="*100)
    out("Similar = same group: %s.   Change is read as %s -> %s" % (" ".join(GROUPS), LBL_REF, LBL_MOB))
    out("%-4s %-9s %-9s %-14s %s" % ("Pos", LBL_REF, LBL_MOB, "Category", "Change"))
    for p in pos:
        out("%-4s %-9s %-9s %-14s %s" % ("P%d" % p["n"], lab(p["a"]), lab(p["b"]), p["cat"], change(p["x"], p["y"])))
    cnt = {c: sum(p["cat"] == c for p in pos) for c in ("IDENTICAL","SIMILAR","DIFFERENT","NO EQUIVALENT")}
    out("Totals: " + ", ".join("%s %d" % (k, v) for k, v in cnt.items()))

    def pr(aa):
        if aa is None: return ("-",)*6
        return ("%.1f" % KD[aa], pclass(aa), "Y" if aa in AROM else "-", ionis(aa),
                str(DON.get(aa,0)), str(ACC.get(aa,0)))
    out("\n" + "="*100 + "\nTABLE 4. PROPERTIES OF THE RESIDUES AT EACH POSITION (side chains only)\n" + "="*100)
    out("KD = Kyte-Doolittle hydrophobicity (+ve hydrophobic).  Ar = aromatic.  Ion = ionisable.  D/A = side-chain H-bond donors/acceptors.")
    h = "%-4s | %-8s %5s %-12s %-3s %-10s %-2s %-2s | %-8s %5s %-12s %-3s %-10s %-2s %-2s"
    out(h % ("Pos", LBL_REF, "KD", "Class", "Ar", "Ion", "D", "A", LBL_MOB, "KD", "Class", "Ar", "Ion", "D", "A"))
    for p in pos:
        out(h % (("P%d" % p["n"], lab(p["a"])) + pr(p["x"]) + (lab(p["b"]),) + pr(p["y"])))
    out("Backbone: every residue except Pro also has a backbone N-H donor, and every residue a C=O acceptor (not counted above).")

    out("\n" + "="*100 + "\nTABLE 5. POCKET TOTALS\n" + "="*100)
    sets = [("%s positions" % LBL_REF, totals([p["x"] for p in pos])),
            ("%s positions" % LBL_MOB, totals([p["y"] for p in pos])),
            ("%s pocket only" % LBL_REF, totals([res[REF][k]["aa"] for k in pk[REF]])),
            ("%s pocket only" % LBL_MOB, totals([res[MOB][k]["aa"] for k in pk[MOB]]))]
    out("%-32s" % "property" + "".join("%-18s" % s[0] for s in sets))
    for k in sets[0][1]:
        out("%-32s" % k + "".join("%-18s" % str(s[1].get(k, "")) for s in sets))
    out("('positions' columns compare like with like, using the residues at every matched position.)")

    # ---- 5. files
    txt = os.path.join(OUTDIR, "pocket_report.txt")
    with open(txt, "w", encoding="utf-8") as fh: fh.write("\n".join(L))
    csvf = os.path.join(OUTDIR, "pocket_positions.csv")
    with open(csvf, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        hd = ["pos", LBL_REF, LBL_MOB, "in_pocket_of", "CA_CA_A", "sidechain_centroid_A", "paired_by", "category", "change"]
        for l in (LBL_REF, LBL_MOB): hd += ["%s_%s" % (l, c) for c in ("KD","class","aromatic","ionisable","donors","acceptors")]
        w.writerow(hd)
        for p in pos:
            w.writerow(["P%d" % p["n"], lab(p["a"]), lab(p["b"]), p["where"],
                        "" if p["dca"] is None else "%.2f" % p["dca"], "" if p["dsc"] is None else "%.2f" % p["dsc"],
                        p["how"], p["cat"], change(p["x"], p["y"])] + list(pr(p["x"])) + list(pr(p["y"])))
    print("\nWritten:\n  %s\n  %s" % (txt, csvf))

    # ---- 6. display
    col = {"IDENTICAL":"green", "SIMILAR":"yellow", "DIFFERENT":"red", "NO EQUIVALENT":"gray50"}
    cmd.hide("everything", "all")
    cmd.show("cartoon", "%s or %s" % (REF, MOB)); cmd.set("cartoon_transparency", 0.7)
    cmd.color("lightblue", REF); cmd.color("wheat", MOB)
    cmd.show("sticks", "%s and (resn ADP or resn SO4)" % REF); cmd.color("magenta", "%s and resn ADP and elem C" % REF)
    cmd.set("label_size", 14)
    for p in pos:
        for rec, tag in ((p["a"], LBL_REF), (p["b"], LBL_MOB)):
            if not rec: continue
            s = "byres (%s and index %d)" % (rec["obj"], rec["idx"])
            cmd.show("sticks", s); cmd.color(col[p["cat"]], s + " and elem C")
            cmd.label("%s and index %d" % (rec["obj"], rec["idx"]), repr("P%d %s" % (p["n"], lab(rec))))
        if DRAW_LINES and p["a"] and p["b"]:
            cmd.distance("pair_lines", "%s and index %d" % (REF, p["a"]["idx"]), "%s and index %d" % (MOB, p["b"]["idx"]))
    if DRAW_LINES:
        cmd.hide("labels", "pair_lines"); cmd.color("gray60", "pair_lines")
    cmd.zoom("probe", 10); cmd.deselect()
    print("Colours: green=identical, yellow=similar, red=different, grey=no equivalent. Labels: P<position> residue.")

main()