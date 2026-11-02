#!/usr/bin/env python3
"""Definitive, internally-consistent study-bias-matched null.
All sets restricted to the annotated background pool; one deterministic pass."""
import sys, os, random, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pyhpo import Ontology

GS_DIR, FIG_DIR, RES_DIR = sys.argv[1], sys.argv[2], sys.argv[3]
N_PERM, SEED = 10000, 42

def load(name):
    s=set()
    with open(os.path.join(GS_DIR,name+".txt")) as fh:
        for l in fh:
            g=l.strip()
            if g and g!="-": s.add(g)
    return s

Ontology()
from pyhpo import Ontology as O
ann={g.name:len(g.hpo) for g in O.genes}

background=load("background_genes")
pool=sorted(g for g in background if g in ann)          # deterministic order
pool_set=set(pool)
bg_counts=np.array([ann[g] for g in pool])

# restrict every set to the pool
macro=load("macrocephaly") & pool_set
comps={c: load(c) & pool_set for c in ["microcephaly","autism","seizure"]}

# bias bins (deciles of pool annotation load)
edges=np.quantile(bg_counts,np.linspace(0,1,11)); edges[-1]+=1
def bin_of(c): return int(np.clip(np.searchsorted(edges,c,side="right")-1,0,9))
bins={i:[] for i in range(10)}
for g in pool: bins[bin_of(ann[g])].append(g)
mbc={i:0 for i in range(10)}
for g in macro: mbc[bin_of(ann[g])]+=1

def run_test(target_set, rng):
    """target_set fixed; resample macro-matched & size-only; return stats."""
    obs=len(macro & target_set)
    nm=np.empty(N_PERM); ns=np.empty(N_PERM)
    n=len(macro)
    for i in range(N_PERM):
        pm=set()
        for b,k in mbc.items():
            if k: pm.update(rng.sample(bins[b],k))
        nm[i]=len(pm & target_set)
        ns[i]=len(set(rng.sample(pool,n)) & target_set)
    p_m=(np.sum(nm>=obs)+1)/(N_PERM+1)
    p_s=(np.sum(ns>=obs)+1)/(N_PERM+1)
    return dict(obs=obs, s_mean=ns.mean(), s_fold=obs/ns.mean(), s_p=p_s,
                m_mean=nm.mean(), m_fold=obs/nm.mean(), m_p=p_m)

rng=random.Random(SEED)
rows=[]
for c in ["microcephaly","autism","seizure"]:
    r=run_test(comps[c], rng)
    rows.append((f"macrocephaly vs {c}", r["obs"],
                 round(r["s_mean"],1), round(r["s_fold"],2), r["s_p"],
                 round(r["m_mean"],1), round(r["m_fold"],2), r["m_p"]))
    print(f"vs {c:12s}: obs={r['obs']:4d} | size-only fold={r['s_fold']:.2f} p={r['s_p']:.3g}"
          f" | bias-matched fold={r['m_fold']:.2f} p={r['m_p']:.3g}")

# four-way shared core
other3=comps["microcephaly"] & comps["autism"] & comps["seizure"]
shared=macro & other3
r4=run_test(other3, rng)
print(f"\nfour-way core: obs={r4['obs']} genes | size-only fold={r4['s_fold']:.2f} p={r4['s_p']:.3g}"
      f" | bias-matched fold={r4['m_fold']:.2f} p={r4['m_p']:.3g}")
sh_counts=np.array([ann[g] for g in shared])
print(f"shared-core study load: median {np.median(sh_counts):.0f} HPO terms vs background median {np.median(bg_counts):.0f}")

with open(os.path.join(RES_DIR,"overlap_null_biasmatched.csv"),"w",newline="") as fh:
    w=csv.writer(fh)
    w.writerow(["comparison","observed_overlap","sizeonly_null_mean","sizeonly_fold","sizeonly_p",
                "biasmatched_null_mean","biasmatched_fold","biasmatched_p"])
    w.writerows(rows)
    w.writerow([f"four-way shared core (n={r4['obs']})", r4["obs"],
                round(r4["s_mean"],1),round(r4["s_fold"],2),r4["s_p"],
                round(r4["m_mean"],1),round(r4["m_fold"],2),r4["m_p"]])

with open(os.path.join(RES_DIR,"shared_core_genes.txt"),"w") as fh:
    fh.write(f"Macrocephaly genes shared with microcephaly, autism AND seizure (n={len(shared)}):\n")
    fh.write(", ".join(sorted(shared))+"\n")
    fh.write(f"\nMedian HPO annotations: shared core={np.median(sh_counts):.0f}, background={np.median(bg_counts):.0f}\n")

# ---- figure ----
labels=["microcephaly","autism","seizure"]
rawf=[rw[3] for rw in rows]; bmf=[rw[6] for rw in rows]; bmp=[rw[7] for rw in rows]
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(11,4.5))
x=np.arange(3); w=0.36
ax1.bar(x-w/2,rawf,w,label="Size-only null\n(≈ hypergeometric)",color="#9ecae1")
ax1.bar(x+w/2,bmf,w,label="Study-bias-matched null",color="#d95f4c")
ax1.axhline(1,ls="--",lw=1,color="#555")
ax1.set_xticks(x); ax1.set_xticklabels([f"vs {c}" for c in labels])
ax1.set_ylabel("Fold enrichment of overlap with macrocephaly")
ax1.set_title("Most overlap is study bias, not biology")
for i,p in enumerate(bmp):
    tag="p<0.001" if p<1e-3 else f"p={p:.2f}"
    ax1.text(x[i]+w/2,bmf[i]+0.04,tag,ha="center",fontsize=9)
ax1.legend(fontsize=8,frameon=False); ax1.set_ylim(0,max(rawf)*1.25)
ax2.boxplot([bg_counts,sh_counts],tick_labels=["All HPO\ngenes",f"{len(shared)} shared\ncore genes"],
            showfliers=False,widths=0.5,medianprops=dict(color="#d95f4c",lw=2))
ax2.set_ylabel("HPO terms annotated per gene (study load)")
ax2.set_title("Shared-core genes are far more heavily studied")
ax2.set_ylim(0,np.quantile(sh_counts,0.95)*1.1)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR,"bias_matched_overlap.png"),dpi=150)
print("\nsaved figure + results (consistent single pass)")
