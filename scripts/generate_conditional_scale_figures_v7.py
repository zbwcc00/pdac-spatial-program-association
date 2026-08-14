"""Generate post-hoc sensitivity figures and an all-patient display atlas."""
from __future__ import annotations

import csv
import importlib.util
import json
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
from PIL import Image

PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "figures/submission_v6"
SUPP = PROJECT / "manuscript/supplementary_v3_final_mask_moran"
DATA = PROJECT / "data/03_results/conditional_scale_joint_sensitivities_v7"
V2 = PROJECT / "data/03_results/unified_primary_pipeline_v2"
SPATIAL = PROJECT / "data/01_unpacked/spatial/GSE278687_spatial"
BLUE, ORANGE, TEAL, PURPLE, GREY = "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#6B7280"
plt.rcParams.update({"font.family":"Arial", "font.size":8, "axes.linewidth":.7, "pdf.fonttype":42, "ps.fonttype":42})


def read_tsv(path):
    with path.open(encoding="utf-8") as handle: return list(csv.DictReader(handle, delimiter="\t"))


def save(fig, stem):
    OUT.mkdir(parents=True, exist_ok=True); SUPP.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(); fig.savefig(OUT / f"{stem}.png", dpi=400, bbox_inches="tight"); fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight"); plt.close(fig)
    for suffix in (".png", ".pdf"): shutil.copy2(OUT / f"{stem}{suffix}", SUPP / f"{stem}{suffix}")


def panel_sensitivity():
    rows = read_tsv(DATA / "cohort_conditional_scale_summary.tsv")
    masks = ["whole_tissue", "DC_core_above_median", "DC_core_above_mean", "DC_core_above_q75"]
    labels = ["Whole\ntissue", "DC-core\n> median", "DC-core\n> mean", "DC-core\n> Q75"]
    fig, axes = plt.subplots(2, 2, figsize=(7.3,5.5), sharey="row")
    for col, cohort in enumerate(("GSE278687","GSE277116")):
        for row_idx, adjustment in enumerate(("without_DC_core_adjustment","with_DC_core_adjustment")):
            ax=axes[row_idx,col]
            subset=[r for r in rows if r["cohort"]==cohort and r["k"]=="6" and r["mregDC_residualization"]==adjustment]
            ordered=[next(r for r in subset if r["mask_definition"]==m) for m in masks]
            y=np.asarray([float(r["median_local_r"]) for r in ordered]); lo=np.asarray([float(r["bootstrap_ci_low"]) for r in ordered]); hi=np.asarray([float(r["bootstrap_ci_high"]) for r in ordered])
            ax.errorbar(range(4), y, yerr=[y-lo,hi-y], marker="o", color=BLUE if cohort=="GSE278687" else ORANGE, lw=1.5,capsize=2)
            ax.set_xticks(range(4),labels); ax.axhline(0,color=GREY,lw=.7); ax.grid(axis="y",alpha=.18); ax.spines[["top","right"]].set_visible(False)
            ax.set_title(f"{cohort}: {'without' if row_idx==0 else 'with'} DC-core residualization",loc="left",fontsize=8.7,weight="bold")
            if col==0: ax.set_ylabel("Median local r (95% bootstrap CI)")
    fig.suptitle("Post-hoc conditioning and DC-core-mask threshold sensitivity (k=6)",y=1.02,fontsize=10.6,weight="bold")
    save(fig,"Supplementary_Figure_S7_conditioning_and_mask_sensitivity")


def panel_scale_vif():
    rows=read_tsv(DATA / "cohort_conditional_scale_summary.tsv")
    vif=read_tsv(DATA / "joint_model_vif_summary.tsv")
    fig, axes=plt.subplots(1,2,figsize=(7.3,3.15))
    for cohort,color in (("GSE278687",BLUE),("GSE277116",ORANGE)):
        subset=[r for r in rows if r["cohort"]==cohort and r["mregDC_residualization"]=="with_DC_core_adjustment" and r["mask_definition"]=="DC_core_above_mean"]
        subset=sorted(subset,key=lambda x:int(x["k"])); x=np.arange(len(subset)); y=np.array([float(r["median_local_r"]) for r in subset]); lo=np.array([float(r["bootstrap_ci_low"]) for r in subset]); hi=np.array([float(r["bootstrap_ci_high"]) for r in subset])
        axes[0].errorbar(x,y,yerr=[y-lo,hi-y],marker="o",color=color,lw=1.5,capsize=2,label=cohort)
    axes[0].set_xticks(range(3),["k=4","k=6 (primary)","k=12"]); axes[0].set_ylabel("Median local r (95% bootstrap CI)"); axes[0].axhline(0,color=GREY,lw=.7); axes[0].grid(axis="y",alpha=.18); axes[0].spines[["top","right"]].set_visible(False); axes[0].legend(frameon=False,fontsize=7)
    names=["Tfh-like","Broad T","Non-Tfh CD4","Treg","Exhausted CD8","cDC","Macrophage"]
    lookup={"Tfh-like":"Tfh_like","Broad T":"broad_T_nonoverlap","Non-Tfh CD4":"non_Tfh_CD4_nonoverlap","Treg":"Treg_nonoverlap","Exhausted CD8":"exhausted_CD8_nonoverlap","cDC":"cDC_nonoverlap","Macrophage":"macrophage_nonoverlap"}
    x=np.arange(len(names)); width=.36
    for off,cohort,color in ((-width/2,"GSE278687",BLUE),(width/2,"GSE277116",ORANGE)):
        values=[float(next(r for r in vif if r["cohort"]==cohort and r["predictor"]==lookup[n])["median_VIF"]) for n in names]
        axes[1].bar(x+off,values,width,color=color,label=cohort)
    axes[1].axhline(5,color=GREY,ls="--",lw=.8,label="VIF=5 reference"); axes[1].set_xticks(x,names,rotation=45,ha="right"); axes[1].set_ylabel("Patient/sample median VIF"); axes[1].spines[["top","right"]].set_visible(False); axes[1].legend(frameon=False,fontsize=6.8)
    axes[0].set_title("Neighbourhood-size sensitivity",loc="left",weight="bold"); axes[1].set_title("Joint model collinearity audit",loc="left",weight="bold")
    fig.suptitle("Sensitivity to local scale and descriptive joint-model diagnostics",y=1.02,fontsize=10.6,weight="bold")
    save(fig,"Supplementary_Figure_S8_neighbour_scale_and_joint_model_diagnostics")


def load_primary():
    spec=importlib.util.spec_from_file_location("v2_fig_atlas",PROJECT / "scripts/run_unified_primary_pipeline_v2.py"); mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    mod.PROJECT=PROJECT; mod.FROZEN=json.loads((PROJECT / "data/03_results/GSE154778_program_freeze/frozen_programs.json").read_text(encoding="utf-8")); mod.G278=PROJECT / "data/01_unpacked/spatial/GSE278687"; mod.G278_COORDS=SPATIAL
    return mod


def patient_atlas():
    mod=load_primary(); patient_rows=read_tsv(V2 / "GSE278687_v2_per_patient.tsv"); section_rows=read_tsv(V2 / "GSE278687_v2_per_section.tsv")
    patient_effect={r["patient"]:float(r["primary_mregDC_Tfh_local_r"]) for r in patient_rows}
    # One section per patient: minimum absolute distance to that patient's
    # patient-level median, so repeated sections do not receive visual weight.
    selected=[]
    for patient,effect in patient_effect.items():
        candidates=[r for r in section_rows if r["unit"]==patient]
        chosen=min(candidates,key=lambda r:abs(float(r["primary_mregDC_Tfh_local_r"])-effect))
        selected.append((patient,chosen["sample"],effect))
    selected.sort(key=lambda x:x[2])
    prepared=[]
    for patient,sample,effect in selected:
        _,matrix,genes,coords,_,_,_=mod.read_g278(PROJECT / "data/01_unpacked/spatial/GSE278687" / f"{sample}_filtered_feature_bc_matrix.h5")
        scores=mod.add_mreg_no_ccl19(mod.make_scores(matrix,genes),matrix,genes); mreg,tfh,mask,neigh,_=mod.prepared_fields(scores,coords)
        joint=np.minimum(mod.z(mod.local_field(mod.z(mreg),neigh)),mod.z(mod.local_field(mod.z(tfh),neigh)))
        spatial=SPATIAL/sample/"spatial"; scale=json.loads((spatial/"scalefactors_json.json").read_text())["tissue_hires_scalef"]; image=np.asarray(Image.open(spatial/"tissue_hires_image.png").convert("RGB")); prepared.append((patient,sample,effect,image,coords*scale,joint,mask))
    limit=float(np.quantile(np.abs(np.concatenate([p[5][p[6]] for p in prepared])),.98))
    fig,axes=plt.subplots(3,6,figsize=(9.2,4.7)); axes=axes.ravel()
    audit=[]
    scatter = None
    for ax,item in zip(axes,prepared):
        patient,sample,effect,image,coords,joint,mask=item; ax.imshow(image,origin="upper"); scatter=ax.scatter(coords[mask,0],coords[mask,1],c=joint[mask],s=2.2,cmap="viridis",norm=TwoSlopeNorm(vcenter=0,vmin=-limit,vmax=limit),linewidths=0,alpha=.9); ax.set_title(f"{patient}  r={effect:.3f}",fontsize=7,weight="bold",pad=2); ax.axis("off"); audit.append({"patient":patient,"representative_section":sample,"patient_level_primary_local_r":effect,"selection_rule":"section nearest patient-level median; all 18 patients shown"})
    fig.subplots_adjust(left=.015, right=.90, bottom=.04, top=.89, wspace=.08, hspace=.16)
    colorbar=fig.colorbar(scatter,cax=fig.add_axes([.92,.18,.012,.60])); colorbar.set_label("Joint residual min(z)",fontsize=7); colorbar.ax.tick_params(labelsize=6)
    fig.suptitle("GSE278687 all-patient atlas: DC-core-mask local joint residual min(z) field",y=.97,fontsize=10,weight="bold"); fig.savefig(OUT/"Supplementary_Figure_S6_GSE278687_all_patient_joint_field_atlas.png",dpi=400,bbox_inches="tight"); fig.savefig(OUT/"Supplementary_Figure_S6_GSE278687_all_patient_joint_field_atlas.pdf",bbox_inches="tight"); plt.close(fig)
    for suffix in (".png",".pdf"): shutil.copy2(OUT/f"Supplementary_Figure_S6_GSE278687_all_patient_joint_field_atlas{suffix}",SUPP/f"Supplementary_Figure_S6_GSE278687_all_patient_joint_field_atlas{suffix}")
    with (OUT/"Supplementary_Figure_S6_selection_audit.tsv").open("w",newline="",encoding="utf-8") as f: w=csv.DictWriter(f,fieldnames=list(audit[0]),delimiter="\t"); w.writeheader();w.writerows(audit)
    shutil.copy2(OUT/"Supplementary_Figure_S6_selection_audit.tsv",SUPP/"Supplementary_Figure_S6_selection_audit.tsv")


if __name__=="__main__":
    panel_sensitivity(); panel_scale_vif(); patient_atlas(); print(OUT)
