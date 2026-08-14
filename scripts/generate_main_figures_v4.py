"""Generate v4 figures with complementary global- and mask-Moran sensitivities."""
from __future__ import annotations
import csv
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

PROJECT = Path(__file__).resolve().parents[1]
LOCKED = PROJECT / "data/03_results/unified_primary_pipeline_v2"
V3 = PROJECT / "data/03_results/spatial_nulls_v3_999"
V4 = PROJECT / "data/03_results/spatial_nulls_v4_mask_moran"
OUT = PROJECT / "figures/main_v4"
OUT.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, TEAL, PURPLE, GREY = "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#6B7280"
plt.rcParams.update({"font.family":"Arial","font.size":8,"axes.linewidth":.7,"pdf.fonttype":42,"ps.fonttype":42})

def read(path):
    with path.open(encoding="utf-8") as handle: return list(csv.DictReader(handle, delimiter="\t"))
def num(row,key): return float(row[key])
def save(fig,name):
    fig.tight_layout()
    for ext in ("png","pdf"): fig.savefig(OUT / f"{name}.{ext}",dpi=400,bbox_inches="tight")
    plt.close(fig)
def add_null(ax, specs, xlabel):
    for i,(label,row,color,note) in enumerate(specs):
        ax.plot([num(row,"null_ci_low"),num(row,"null_ci_high")],[i,i],color=color,lw=7,alpha=.45,solid_capstyle="round")
        observed=num(row,"observed"); ax.scatter(observed,i,color=BLUE,edgecolor="white",linewidth=.6,s=50,zorder=3)
        ax.text(observed+.010,i+.07,note,fontsize=6.25,va="bottom")
    ax.axvline(0,color=GREY,lw=.8); ax.set_yticks(range(len(specs)),[item[0] for item in specs]); ax.set_xlabel(xlabel)
    ax.grid(axis="x",alpha=.2); ax.spines[["top","right"]].set_visible(False)

def fig1():
    patients=sorted(read(LOCKED / "GSE278687_v2_per_patient.tsv"),key=lambda row:num(row,"primary_mregDC_Tfh_local_r"))
    block=read(V3 / "GSE278687_v3_patient_block_null_999.tsv")[0]
    global_row={row["cohort"]:row for row in read(V3 / "v3_combined_1000_graph_nulls.tsv")}["GSE278687"]
    mask=read(V4 / "GSE278687_v4_patient_mask_moran_combined_1000.tsv")[0]
    fig,axes=plt.subplots(1,2,figsize=(7.35,3.5),gridspec_kw={"width_ratios":[1.25,1.15]})
    ax=axes[0]; y=np.arange(len(patients)); values=[num(row,"primary_mregDC_Tfh_local_r") for row in patients]
    ax.scatter(values,y,color=BLUE,edgecolor="white",linewidth=.5,s=37,zorder=3); ax.axvline(0,color=GREY,lw=.8); ax.axvline(np.median(values),color=ORANGE,lw=1.3,ls="--")
    ax.set_yticks(y,[row["patient"] for row in patients]); ax.set_xlabel("Local program correlation (r)"); ax.set_title("GSE278687: patient-level primary effect",loc="left",weight="bold")
    ax.text(.98,.02,"Median r = 0.418\n17/18 positive",transform=ax.transAxes,ha="right",va="bottom",fontsize=7.2,bbox={"facecolor":"white","edgecolor":"none","alpha":.9,"pad":1.5}); ax.grid(axis="x",alpha=.2); ax.spines[["top","right"]].set_visible(False)
    specs=[("Array-block null (primary)",block,ORANGE,"999 draws; 0/999\nadd-one P = 0.001"),("Global-Moran graph sensitivity",global_row,TEAL,"1,000 draws; 0/1,000\nadd-one Monte Carlo P = 0.000999"),("Mask-Moran graph sensitivity",mask,PURPLE,"1,000 draws; 0/1,000\nadd-one Monte Carlo P = 0.000999")]
    add_null(axes[1],specs,"Cohort statistic: median patient r"); axes[1].set_title("Null models with stated calibration targets",loc="left",weight="bold"); axes[1].set_ylim(-.6,2.65)
    save(fig,"Figure1_v4_primary_and_spatial_nulls")

def fig2():
    summary278={row["endpoint"]:row for row in read(LOCKED / "GSE278687_v2_summary.tsv")}; summary277={row["endpoint"]:row for row in read(LOCKED / "GSE277116_v2_summary.tsv") if row["analysis_unit"]=="sample (main 18)"}
    endpoints=[("primary_mregDC_Tfh_local_r","Primary local r"),("mregDC_Tfh_after_nonoverlap_broad_T_adjustment_r","After non-overlap\nbroad-T adjustment"),("mregDC_no_CCL19_Tfh_local_r","mregDC score\nwithout CCL19"),("joint_model_Tfh_standardized_beta","Descriptive joint-model\nTfh standardized beta")]
    fig,axes=plt.subplots(1,2,figsize=(7.1,3.35),gridspec_kw={"width_ratios":[1.18,1]}); ax=axes[0]; ypos=np.arange(len(endpoints))
    for offset,summary,label,color in [(-.15,summary278,"GSE278687 patients",BLUE),(.15,summary277,"GSE277116 samples",ORANGE)]:
        med=[num(summary[key],"median") for key,_ in endpoints]; lo=[num(summary[key],"bootstrap_ci_low") for key,_ in endpoints]; hi=[num(summary[key],"bootstrap_ci_high") for key,_ in endpoints]
        ax.errorbar(med,ypos+offset,xerr=[np.subtract(med,lo),np.subtract(hi,med)],fmt="o",color=color,capsize=2,label=label,ms=4)
    ax.axvline(0,color=GREY,lw=.8); ax.set_yticks(ypos,[label for _,label in endpoints]); ax.invert_yaxis(); ax.set_xlabel("Effect estimate (95% bootstrap CI)"); ax.set_title("Robustness of local association",loc="left",weight="bold"); ax.legend(frameon=False,loc="lower right",fontsize=7); ax.grid(axis="x",alpha=.2); ax.spines[["top","right"]].set_visible(False)
    ax=axes[1]; comps=[("broad_T_nonoverlap","Broad T"),("non_Tfh_CD4_nonoverlap","Non-Tfh CD4"),("Treg_nonoverlap","Treg"),("exhausted_CD8_nonoverlap","Exhausted CD8"),("cDC_nonoverlap","Conventional DC"),("macrophage_nonoverlap","Macrophage")]; ypos=np.arange(len(comps))
    for offset,summary,color in [(-.15,summary278,BLUE),(.15,summary277,ORANGE)]:
        keys=[f"mregDC_{key}_competitive_r" for key,_ in comps]; med=[num(summary[key],"median") for key in keys]; lo=[num(summary[key],"bootstrap_ci_low") for key in keys]; hi=[num(summary[key],"bootstrap_ci_high") for key in keys]
        ax.errorbar(med,ypos+offset,xerr=[np.subtract(med,lo),np.subtract(hi,med)],fmt="o",color=color,capsize=2,ms=4)
    ax.axvline(0,color=GREY,lw=.8); ax.set_yticks(ypos,[label for _,label in comps]); ax.invert_yaxis(); ax.set_xlabel("Reference-program local correlation (95% bootstrap CI)"); ax.set_title("Non-overlapping reference programs",loc="left",weight="bold"); ax.grid(axis="x",alpha=.2); ax.spines[["top","right"]].set_visible(False)
    save(fig,"Figure2_v4_robustness_and_competition")

def fig3():
    packages=[row for row in read(LOCKED / "GSE277116_v2_per_package.tsv") if row["analysis_set"]=="main_18"]; packages.sort(key=lambda row:num(row,"primary_mregDC_Tfh_local_r"))
    global_row={row["cohort"]:row for row in read(V3 / "v3_combined_1000_graph_nulls.tsv")}["GSE277116"]; mask=read(V4 / "GSE277116_v4_sample_mask_moran_combined_1000.tsv")[0]
    fig,axes=plt.subplots(1,2,figsize=(7.35,3.5),gridspec_kw={"width_ratios":[1.3,1.1]}); ax=axes[0]; y=np.arange(len(packages)); colors=[ORANGE if row["stratum"]=="ffpe" else BLUE for row in packages]; values=[num(row,"primary_mregDC_Tfh_local_r") for row in packages]
    ax.scatter(values,y,color=colors,edgecolor="white",lw=.5,s=37,zorder=3); ax.axvline(0,color=GREY,lw=.8); ax.axvline(np.median(values),color=TEAL,lw=1.2,ls="--"); ax.set_yticks(y,[row["sample"] for row in packages]); ax.set_xlabel("Local program correlation (r)"); ax.set_title("GSE277116: sample-level technical replication",loc="left",weight="bold"); ax.text(.02,.02,"Median r = 0.414\n18/18 positive",transform=ax.transAxes,fontsize=7.2,bbox={"facecolor":"white","edgecolor":"none","alpha":.9,"pad":1.5}); ax.grid(axis="x",alpha=.2); ax.spines[["top","right"]].set_visible(False)
    specs=[("Global-Moran graph sensitivity",global_row,TEAL,"0/1,000\nadd-one Monte Carlo P = 0.000999"),("Mask-Moran graph sensitivity",mask,PURPLE,"0/1,000\nadd-one Monte Carlo P = 0.000999")]
    add_null(axes[1],specs,"Cohort statistic: median sample r"); axes[1].set_title("Complementary spatial sensitivities",loc="left",weight="bold"); axes[1].set_ylim(-.55,1.55)
    save(fig,"Figure3_v4_external_technical_replication")

if __name__=="__main__": fig1(); fig2(); fig3()
