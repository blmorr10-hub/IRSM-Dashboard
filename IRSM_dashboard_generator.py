"""
Immune Risk Scoring Module — Dashboard Generator
-------------------------------------------------
Biomedical Product Design — Team 13
Run this script to:
  1. Compute the full IRSM score from biomarker + context inputs
  2. Automatically generate a styled HTML dashboard with the real output values
  3. Open the dashboard in your browser

Usage:
  python3 IRSM_dashboard_generator.py

To change the patient, edit the BIOMARKERS and CONTEXT dictionaries below.
"""

import webbrowser
import os
from datetime import datetime


# ── STEP 1: BIOMARKER AND CONTEXT INPUTS ─────────────────────────────────────
# Edit these values to model a different patient profile.

BIOMARKERS = {
    # Innate Inflammation (pg/mL)
    "IL6":       22.0,
    "CRP":       14.0,
    "TNFa":      10.0,
    # Complement (ng/mL)
    "C3a":       420.0,
    "C5b9":      310.0,
    # T-Cell Activation (pg/mL)
    "IFNg":      14.0,
    "GranzymeB": 280.0,
    "IL17":      30.0,
    # Antibody / Fc
    "DSA":       1200.0,   # MFI
    "IgG":       1800.0,   # mg/dL
    # Fibrosis Remodeling (pg/mL or ng/mL)
    "TGFb":      1400.0,
    "MMP2":      240.0,
    "CTGF":      180.0,
}

CONTEXT = {
    "hla_mismatch":        3,
    "implant_age_months":  18,
    "prior_rejection":     False,
    "immunosuppression":   "low",      # "high" | "moderate" | "low" | "none"
    "biopsy_grade":        "borderline",  # "none"|"borderline"|"1A"|"1B"|"2"|"3"
    "implant_type":        "Titanium alloy",
}


# ── STEP 2: REFERENCE RANGES AND WEIGHTS ─────────────────────────────────────

REFERENCE_RANGES = {
    "IL6":       {"normal": 7,    "critical": 30},
    "CRP":       {"normal": 5,    "critical": 20},
    "TNFa":      {"normal": 8.1,  "critical": 25},
    "C3a":       {"normal": 200,  "critical": 700},
    "C5b9":      {"normal": 250,  "critical": 900},
    "IFNg":      {"normal": 5,    "critical": 25},
    "GranzymeB": {"normal": 100,  "critical": 500},
    "IL17":      {"normal": 20,   "critical": 100},
    "DSA":       {"normal": 500,  "critical": 2000},
    "IgG":       {"normal": 1500, "critical": 3000},
    "TGFb":      {"normal": 1000, "critical": 4000},
    "MMP2":      {"normal": 200,  "critical": 700},
    "CTGF":      {"normal": 150,  "critical": 600},
}

# Intra-pathway biomarker weights (wₖⱼ) — must sum to 1 per pathway
PATHWAY_MARKER_WEIGHTS = {
    "Innate Inflammation": {"IL6": 0.40, "CRP": 0.30, "TNFa": 0.30},
    "Complement":          {"C3a": 0.50, "C5b9": 0.50},
    "T-Cell Activation":   {"IFNg": 0.35, "GranzymeB": 0.40, "IL17": 0.25},
    "Antibody/Fc":         {"DSA": 0.60, "IgG": 0.40},
    "Fibrosis Remodeling": {"TGFb": 0.40, "MMP2": 0.35, "CTGF": 0.25},
}

# Pathway-to-process weights (vᵢₖ) — must sum to 1 per process
PROCESS_PATHWAY_WEIGHTS = {
    "Inflammation": {"Innate Inflammation": 0.60, "Complement": 0.40},
    "Tissue Injury": {"T-Cell Activation": 0.55, "Antibody/Fc": 0.45},
    "Fibrosis":      {"Fibrosis Remodeling": 1.00},
}

# Process-level weights (Wᵢ) — must sum to 1
PROCESS_WEIGHTS = {"Inflammation": 0.40, "Tissue Injury": 0.40, "Fibrosis": 0.20}

# Context modifier coefficients (cₗ · mₗ)
CONTEXT_MODIFIERS = {
    "hla_mismatch":        {0: 0.00, 1: 0.02, 2: 0.04, 3: 0.07, 4: 0.10, 5: 0.13, 6: 0.16},
    "implant_age_months":  lambda m: min(m * 0.003, 0.10),
    "prior_rejection":     {True: 0.12, False: 0.00},
    "immunosuppression":   {"high": -0.08, "moderate": 0.00, "low": 0.06, "none": 0.12},
    "biopsy_grade":        {"none": 0.00, "borderline": 0.04, "1A": 0.08, "1B": 0.10, "2": 0.14, "3": 0.18},
}

# Risk tier thresholds (θ₁, θ₂, θ₃)
THRESHOLDS = {"low": 25, "moderate": 50, "mod_high": 70}

ACTIVATION_THRESHOLD = 0.40   # Pₖ ≥ 0.40 → pathway activated


# ── STEP 3: SCORING FUNCTIONS ─────────────────────────────────────────────────

def normalize(marker, value):
    """Normalize a biomarker value to [0, 1] using reference range."""
    ref = REFERENCE_RANGES[marker]
    if value <= ref["normal"]:
        return 0.0
    elif value >= ref["critical"]:
        return 1.0
    return (value - ref["normal"]) / (ref["critical"] - ref["normal"])


def compute_pathway_scores(biomarkers):
    """Tier 1: Pₖ = Σⱼ (wₖⱼ · bₖⱼ)"""
    scores = {}
    for pathway, weights in PATHWAY_MARKER_WEIGHTS.items():
        score = sum(
            weights[m] * normalize(m, biomarkers[m])
            for m in weights if m in biomarkers
        )
        scores[pathway] = {
            "score":     round(score, 4),
            "score_pct": round(score * 100, 1),
            "activated": score >= ACTIVATION_THRESHOLD,
            "markers": {
                m: round(normalize(m, biomarkers[m]) * 100, 1)
                for m in weights if m in biomarkers
            }
        }
    return scores


def compute_process_scores(pathway_scores):
    """Tier 2: Qᵢ = Σₖ∈ᵢ (vᵢₖ · Pₖ)"""
    scores = {}
    for process, weights in PROCESS_PATHWAY_WEIGHTS.items():
        score = sum(
            weights[p] * pathway_scores[p]["score"]
            for p in weights if p in pathway_scores
        )
        scores[process] = round(score, 4)
    return scores


def compute_context_modifier(context):
    """ε = Σₗ (cₗ · mₗ) — additive context adjustment"""
    eps = 0.0
    eps += CONTEXT_MODIFIERS["hla_mismatch"].get(min(context.get("hla_mismatch", 0), 6), 0)
    eps += CONTEXT_MODIFIERS["implant_age_months"](context.get("implant_age_months", 0))
    eps += CONTEXT_MODIFIERS["prior_rejection"].get(context.get("prior_rejection", False), 0)
    eps += CONTEXT_MODIFIERS["immunosuppression"].get(context.get("immunosuppression", "moderate"), 0)
    eps += CONTEXT_MODIFIERS["biopsy_grade"].get(context.get("biopsy_grade", "none"), 0)
    return round(eps, 4)


def compute_rrs(process_scores, context_modifier):
    """
    Tier 3:
      S     = Σᵢ (Wᵢ · Qᵢ)
      S_adj = clip(S + ε, 0, 1)
      RRS   = S_adj × 100
    """
    S     = sum(PROCESS_WEIGHTS[p] * process_scores[p] for p in PROCESS_WEIGHTS)
    S_adj = max(0.0, min(S + context_modifier, 1.0))
    RRS   = round(S_adj * 100, 1)
    return S, S_adj, RRS


def logistic_calibration(S_adj, alpha=6.0, theta=0.50):
    """P(rejection) = 1 / (1 + e^(−α(S_adj − θ)))"""
    import math
    return round(1 / (1 + math.exp(-alpha * (S_adj - theta))), 3)


def classify_risk(rrs):
    if rrs < THRESHOLDS["low"]:      return "Low",          "#00c9a7"
    elif rrs < THRESHOLDS["moderate"]: return "Moderate",   "#45aaf2"
    elif rrs < THRESHOLDS["mod_high"]: return "Moderate–High", "#ff8c42"
    else:                               return "High",       "#ff4f6e"


def suggest_actions(process_scores, pathway_scores, rrs):
    actions = []
    if rrs >= THRESHOLDS["mod_high"]:
        actions.append("HIGH RISK: Recommend urgent biopsy and multidisciplinary review.")
    if process_scores["Inflammation"] * 100 >= 50:
        if pathway_scores["Innate Inflammation"]["activated"]:
            actions.append("Escalate anti-inflammatory protocol — innate pathway dominant.")
        if pathway_scores["Complement"]["activated"]:
            actions.append("Complement pathway activated — evaluate inhibitor therapy.")
    if process_scores["Tissue Injury"] * 100 >= 50:
        if pathway_scores["T-Cell Activation"]["activated"]:
            actions.append("T-cell cytotoxicity elevated — repeat T-cell panel in 4 weeks.")
        if pathway_scores["Antibody/Fc"]["activated"]:
            actions.append("DSA trending — recommend HLA typing review.")
    if process_scores["Fibrosis"] * 100 >= 40:
        actions.append("Early fibrotic remodeling — consider antifibrotic evaluation.")
    else:
        actions.append("Fibrosis score low — no anti-fibrotic adjustment indicated.")
    if not actions:
        actions.append("Continue routine monitoring schedule.")
    return actions


# ── STEP 4: RUN THE PIPELINE ──────────────────────────────────────────────────

pathway_scores    = compute_pathway_scores(BIOMARKERS)
process_scores    = compute_process_scores(pathway_scores)
context_modifier  = compute_context_modifier(CONTEXT)
S, S_adj, RRS     = compute_rrs(process_scores, context_modifier)
p_rejection       = logistic_calibration(S_adj)
risk_label, risk_color = classify_risk(RRS)
actions           = suggest_actions(process_scores, pathway_scores, RRS)

print(f"RRS: {RRS} / 100  |  Category: {risk_label}  |  P(rejection): {p_rejection}")


# ── STEP 5: GENERATE HTML DASHBOARD ──────────────────────────────────────────

def pct(val):
    """Convert 0–1 score to display percentage string."""
    return f"{round(val * 100, 1)}%"

def bar(val, color):
    """Inline progress bar HTML."""
    w = min(int(val * 100), 100)
    return f'<div style="height:6px;background:rgba(255,255,255,0.08);border-radius:99px;overflow:hidden;margin-bottom:8px"><div style="height:100%;width:{w}%;background:{color};border-radius:99px"></div></div>'

def action_items(actions):
    items = ""
    for a in actions:
        items += f'<div style="display:flex;gap:8px;align-items:flex-start;font-size:11px;line-height:1.5;margin-bottom:6px;color:#d4eaf5"><div style="width:5px;height:5px;border-radius:50%;background:#00c9a7;margin-top:5px;flex-shrink:0"></div>{a}</div>'
    return items

# Context modifier breakdown for display
ctx_contributions = {
    "HLA Mismatch":       CONTEXT_MODIFIERS["hla_mismatch"].get(min(CONTEXT.get("hla_mismatch", 0), 6), 0),
    "Implant Age":        CONTEXT_MODIFIERS["implant_age_months"](CONTEXT.get("implant_age_months", 0)),
    "Prior Rejection":    CONTEXT_MODIFIERS["prior_rejection"].get(CONTEXT.get("prior_rejection", False), 0),
    "Immunosuppression":  CONTEXT_MODIFIERS["immunosuppression"].get(CONTEXT.get("immunosuppression", "moderate"), 0),
    "Biopsy Grade":       CONTEXT_MODIFIERS["biopsy_grade"].get(CONTEXT.get("biopsy_grade", "none"), 0),
}

ctx_rows = ""
for label, val in ctx_contributions.items():
    sign = "+" if val >= 0 else ""
    ctx_rows += f"""
    <div style="padding:7px 10px;border-radius:8px;background:rgba(255,255,255,0.03);border:1px solid rgba(255,255,255,0.06)">
      <div style="font-family:'DM Mono',monospace;font-size:9px;color:#5a8099;text-transform:uppercase;letter-spacing:0.12em;margin-bottom:3px">{label}</div>
      <div style="font-size:12px;font-weight:600;color:#d4eaf5">{CONTEXT.get(label.lower().replace(' ','_').replace('-','_'), '—')}</div>
      <div style="font-family:'DM Mono',monospace;font-size:9px;color:#00c9a7">cₗ·mₗ = {sign}{round(val,3)}</div>
    </div>"""

# Pathway rows for Layer 2
pw_colors = {
    "Innate Inflammation": "#ff6b6b",
    "Complement":          "#f7b731",
    "T-Cell Activation":   "#45aaf2",
    "Antibody/Fc":         "#a29bfe",
    "Fibrosis Remodeling": "#fd9644",
}
pw_eq = {
    "Innate Inflammation": "v₁=0.60 → Inflam",
    "Complement":          "v₂=0.40 → Inflam",
    "T-Cell Activation":   "v₃=0.55 → Injury",
    "Antibody/Fc":         "v₄=0.45 → Injury",
    "Fibrosis Remodeling": "v₅=1.00 → Fibrosis",
}

pathway_html = ""
for pw, data in pathway_scores.items():
    col   = pw_colors[pw]
    act   = '<span style="font-size:8px;font-weight:bold;color:#ff4f6e;background:rgba(255,79,110,0.12);border:1px solid rgba(255,79,110,0.3);border-radius:8px;padding:1px 6px;margin-left:6px">ACTIVATED</span>' if data["activated"] else ""
    markers_str = " · ".join(PATHWAY_MARKER_WEIGHTS[pw].keys())
    pathway_html += f"""
    <div style="display:flex;align-items:center;gap:10px;padding:8px 12px;border-radius:8px;background:rgba(255,255,255,0.03);border:1px solid rgba(255,255,255,0.06);margin-bottom:7px">
      <div style="width:8px;height:8px;border-radius:50%;background:{col};flex-shrink:0"></div>
      <div style="flex:1;font-size:12px;font-weight:600;color:#d4eaf5">{pw}{act}</div>
      <div style="font-family:'DM Mono',monospace;font-size:9px;color:#00c9a7">{pw_eq[pw]}</div>
      <div style="font-family:'DM Mono',monospace;font-size:10px;color:#5a8099">{markers_str}</div>
    </div>"""

# Process cards
def process_card(title, eq, bars, sources, bg, border, col):
    bar_html = ""
    for label, val in bars:
        bar_html += f'<div style="font-family:\'DM Mono\',monospace;font-size:9px;color:#5a8099;display:flex;justify-content:space-between;margin-bottom:3px"><span>{label}</span><span>{val}%</span></div>{bar(val/100, col)}'
    src_html = "".join(f'<span style="margin-right:4px;opacity:0.7">↑ {s}</span>' for s in sources)
    return f"""
    <div style="border-radius:14px;padding:18px;background:{bg};border:1px solid {border};color:{col}">
      <div style="font-size:11px;font-weight:700;letter-spacing:0.1em;text-transform:uppercase;font-family:'DM Mono',monospace;margin-bottom:4px">{title}</div>
      <div style="font-family:'DM Mono',monospace;font-size:9px;color:#00c9a7;margin-bottom:10px">{eq}</div>
      {bar_html}
      <div style="font-family:'DM Mono',monospace;font-size:9px;color:#5a8099;margin-top:6px">{src_html}</div>
    </div>"""

inflam_pct  = round(process_scores["Inflammation"] * 100, 1)
injury_pct  = round(process_scores["Tissue Injury"] * 100, 1)
fibrosis_pct = round(process_scores["Fibrosis"] * 100, 1)

innate_pct  = round(pathway_scores["Innate Inflammation"]["score_pct"], 1)
comp_pct    = round(pathway_scores["Complement"]["score_pct"], 1)
tcell_pct   = round(pathway_scores["T-Cell Activation"]["score_pct"], 1)
ab_pct      = round(pathway_scores["Antibody/Fc"]["score_pct"], 1)
fibro_pct   = round(pathway_scores["Fibrosis Remodeling"]["score_pct"], 1)
ecm_pct     = round(fibro_pct * 0.67, 1)   # ECM deposition sub-component

pc_inflam = process_card("Inflammation Score",
    "Q_inflam = v₁·P_innate + v₂·P_complement",
    [("Innate pathway (v₁=0.60)", innate_pct), ("Complement (v₂=0.40)", comp_pct)],
    ["IL-6", "CRP", "C3a"],
    "rgba(255,107,107,0.08)", "rgba(255,107,107,0.3)", "#ff6b6b")

pc_injury = process_card("Tissue Injury Score",
    "Q_injury = v₃·P_tcell + v₄·P_antibody",
    [("T-cell cytotoxicity (v₃=0.55)", tcell_pct), ("Ab-mediated lysis (v₄=0.45)", ab_pct)],
    ["Granzyme B", "DSA"],
    "rgba(247,183,49,0.08)", "rgba(247,183,49,0.3)", "#f7b731")

pc_fibro = process_card("Fibrosis Score",
    "Q_fibrosis = v₅·P_remodel + v₆·P_ecm",
    [("Remodeling (v₅=0.60)", fibro_pct), ("ECM deposition (v₆=0.40)", ecm_pct)],
    ["TGF-β", "CTGF"],
    "rgba(253,150,68,0.08)", "rgba(253,150,68,0.3)", "#fd9644")

# Risk tier rows
tiers = [
    ("Low",          f"RRS < {THRESHOLDS['low']}",                                          "#00c9a7"),
    ("Moderate",     f"{THRESHOLDS['low']} ≤ RRS < {THRESHOLDS['moderate']}",               "#45aaf2"),
    ("Moderate–High",f"{THRESHOLDS['moderate']} ≤ RRS < {THRESHOLDS['mod_high']}",          "#ff8c42"),
    ("High",         f"RRS ≥ {THRESHOLDS['mod_high']}",                                     "#ff4f6e"),
]
tier_html = ""
for t_label, t_range, t_col in tiers:
    is_active = t_label == risk_label
    weight    = "font-weight:600" if is_active else "opacity:0.5"
    arrow     = " ◀" if is_active else ""
    tier_html += f'<div style="display:flex;align-items:center;gap:6px;font-family:\'DM Mono\',monospace;font-size:9px;color:{t_col};{weight};margin-bottom:3px"><div style="width:6px;height:6px;border-radius:50%;background:{t_col};flex-shrink:0"></div>{t_label}: {t_range}{arrow}</div>'

# Rejection process profile bars
profile_bars = [
    ("Inflammation", inflam_pct,  "#ff6b6b"),
    ("T-cell Injury", tcell_pct, "#45aaf2"),
    ("Ab-mediated",  ab_pct,      "#a29bfe"),
    ("Fibrosis",     fibrosis_pct,"#fd9644"),
]
profile_html = ""
for p_label, p_val, p_col in profile_bars:
    profile_html += f"""
    <div style="display:flex;align-items:center;gap:8px;font-size:11px;font-weight:600;margin-bottom:6px">
      <span style="color:{p_col};width:90px">{p_label}</span>
      <div style="height:4px;border-radius:99px;flex:1;background:{p_col};opacity:0.5;width:{p_val}%"></div>
      <div style="font-family:'DM Mono',monospace;font-size:10px;color:#5a8099;width:36px;text-align:right">{p_val}%</div>
    </div>"""

timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

HTML = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>IRSM Dashboard — Generated {timestamp}</title>
<style>
  @import url('https://fonts.googleapis.com/css2?family=DM+Mono:wght@300;400;500&family=Syne:wght@400;600;700;800&display=swap');
  * {{ box-sizing:border-box; margin:0; padding:0; }}
  body {{
    background:#050d12; color:#d4eaf5;
    font-family:'Syne',sans-serif;
    padding:40px 32px 56px; width:1100px;
  }}
  body::before {{
    content:''; position:fixed; inset:0;
    background-image:linear-gradient(rgba(0,201,167,0.03) 1px, transparent 1px),
      linear-gradient(90deg,rgba(0,201,167,0.03) 1px, transparent 1px);
    background-size:40px 40px; pointer-events:none; z-index:0;
  }}
  .wrapper {{ position:relative; z-index:1; max-width:1100px; margin:0 auto; }}
</style>
</head>
<body>
<div class="wrapper">

  <div style="text-align:center;margin-bottom:36px">
    <div style="font-family:'DM Mono',monospace;font-size:11px;letter-spacing:0.2em;color:#00c9a7;text-transform:uppercase;margin-bottom:12px">
      Biomedical Product Design · System 3 Architecture · Generated {timestamp}
    </div>
    <h1 style="font-size:36px;font-weight:800;letter-spacing:-0.02em;background:linear-gradient(135deg,#fff 40%,#0099cc);-webkit-background-clip:text;-webkit-text-fill-color:transparent;background-clip:text;margin-bottom:10px">
      Immune Risk Scoring Module
    </h1>
    <div style="font-family:'DM Mono',monospace;font-size:12px;color:#5a8099;letter-spacing:0.05em">
      Multi-process rejection biology → unified clinical risk score &nbsp;|&nbsp; RRS = S_adj × 100
    </div>
  </div>

  <!-- Pathway tags -->
  <div style="display:flex;justify-content:center;gap:10px;margin-bottom:28px;flex-wrap:wrap">
    {''.join(f'<span style="font-family:\'DM Mono\',monospace;font-size:10px;letter-spacing:0.1em;padding:4px 12px;border-radius:99px;border:1px solid #1a3040;color:#5a8099;text-transform:uppercase">{p}</span>' for p in ["Innate Inflammation","Complement","T-Cell Activation","Antibody / Fc","Fibrosis Remodeling"])}
  </div>

  <!-- LAYER 1 -->
  <div style="font-family:'DM Mono',monospace;font-size:9px;letter-spacing:0.18em;text-transform:uppercase;color:#5a8099;text-align:center;padding:6px 0 10px">
    Layer 1 · Inputs — Normalized biomarker values bₖⱼ ∈ [0,1] + context variables mₗ
  </div>

  <div style="display:grid;grid-template-columns:1fr 1fr;gap:16px;margin-bottom:6px">

    <!-- Biomarker panel -->
    <div style="background:#0a1a24;border:1px solid #1a3040;border-radius:14px;padding:20px 22px;position:relative;overflow:hidden">
      <div style="position:absolute;top:0;left:0;right:0;height:2px;background:linear-gradient(to right,#00c9a7,#0099cc)"></div>
      <div style="font-size:11px;font-weight:700;letter-spacing:0.1em;text-transform:uppercase;color:#00c9a7;margin-bottom:14px;font-family:'DM Mono',monospace">⬡ Pathway-Organized Biomarker Panel</div>
      {pathway_html}
    </div>

    <!-- Context -->
    <div style="background:#0a1a24;border:1px solid #1a3040;border-radius:14px;padding:20px 22px;position:relative;overflow:hidden">
      <div style="position:absolute;top:0;left:0;right:0;height:2px;background:linear-gradient(to right,#0099cc,#a29bfe)"></div>
      <div style="font-size:11px;font-weight:700;letter-spacing:0.1em;text-transform:uppercase;color:#0099cc;margin-bottom:14px;font-family:'DM Mono',monospace">⬡ Patient & Implant Context (mₗ · cₗ) &nbsp; ε = {round(context_modifier,3)}</div>
      <div style="display:grid;grid-template-columns:1fr 1fr;gap:8px">
        <div style="padding:8px 10px;border-radius:8px;background:rgba(255,255,255,0.03);border:1px solid rgba(255,255,255,0.06)">
          <div style="font-family:'DM Mono',monospace;font-size:9px;color:#5a8099;text-transform:uppercase;letter-spacing:0.12em;margin-bottom:3px">Implant Type</div>
          <div style="font-size:12px;font-weight:600;color:#d4eaf5">{CONTEXT['implant_type']}</div>
          <div style="font-family:'DM Mono',monospace;font-size:9px;color:#00c9a7">cₗ·mₗ = +0.000</div>
        </div>
        <div style="padding:8px 10px;border-radius:8px;background:rgba(255,255,255,0.03);border:1px solid rgba(255,255,255,0.06)">
          <div style="font-family:'DM Mono',monospace;font-size:9px;color:#5a8099;text-transform:uppercase;letter-spacing:0.12em;margin-bottom:3px">Duration</div>
          <div style="font-size:12px;font-weight:600;color:#d4eaf5">{CONTEXT['implant_age_months']} months</div>
          <div style="font-family:'DM Mono',monospace;font-size:9px;color:#00c9a7">cₗ·mₗ = +{round(ctx_contributions['Implant Age'],3)}</div>
        </div>
        <div style="padding:8px 10px;border-radius:8px;background:rgba(255,255,255,0.03);border:1px solid rgba(255,255,255,0.06)">
          <div style="font-family:'DM Mono',monospace;font-size:9px;color:#5a8099;text-transform:uppercase;letter-spacing:0.12em;margin-bottom:3px">HLA Mismatch</div>
          <div style="font-size:12px;font-weight:600;color:#d4eaf5">{CONTEXT['hla_mismatch']} antigens</div>
          <div style="font-family:'DM Mono',monospace;font-size:9px;color:#00c9a7">cₗ·mₗ = +{round(ctx_contributions['HLA Mismatch'],3)}</div>
        </div>
        <div style="padding:8px 10px;border-radius:8px;background:rgba(255,255,255,0.03);border:1px solid rgba(255,255,255,0.06)">
          <div style="font-family:'DM Mono',monospace;font-size:9px;color:#5a8099;text-transform:uppercase;letter-spacing:0.12em;margin-bottom:3px">Prior Rejection</div>
          <div style="font-size:12px;font-weight:600;color:#d4eaf5">{'Yes' if CONTEXT['prior_rejection'] else 'None'}</div>
          <div style="font-family:'DM Mono',monospace;font-size:9px;color:#00c9a7">cₗ·mₗ = +{round(ctx_contributions['Prior Rejection'],3)}</div>
        </div>
        <div style="padding:8px 10px;border-radius:8px;background:rgba(255,255,255,0.03);border:1px solid rgba(255,255,255,0.06)">
          <div style="font-family:'DM Mono',monospace;font-size:9px;color:#5a8099;text-transform:uppercase;letter-spacing:0.12em;margin-bottom:3px">Immunosuppression</div>
          <div style="font-size:12px;font-weight:600;color:#d4eaf5">{CONTEXT['immunosuppression'].capitalize()}</div>
          <div style="font-family:'DM Mono',monospace;font-size:9px;color:#00c9a7">cₗ·mₗ = +{round(ctx_contributions['Immunosuppression'],3)}</div>
        </div>
        <div style="padding:8px 10px;border-radius:8px;background:rgba(255,255,255,0.03);border:1px solid rgba(255,255,255,0.06)">
          <div style="font-family:'DM Mono',monospace;font-size:9px;color:#5a8099;text-transform:uppercase;letter-spacing:0.12em;margin-bottom:3px">Biopsy Grade</div>
          <div style="font-size:12px;font-weight:600;color:#d4eaf5">{CONTEXT['biopsy_grade'].capitalize()}</div>
          <div style="font-family:'DM Mono',monospace;font-size:9px;color:#00c9a7">cₗ·mₗ = +{round(ctx_contributions['Biopsy Grade'],3)}</div>
        </div>
      </div>
    </div>
  </div>

  <!-- Equation strip 1 -->
  <div style="width:100%;background:rgba(0,201,167,0.04);border:1px solid rgba(0,201,167,0.15);border-radius:10px;padding:10px 18px;margin:6px 0;font-family:'DM Mono',monospace;font-size:11px;color:#00c9a7;text-align:center">
    <div style="font-size:9px;color:#5a8099;text-transform:uppercase;letter-spacing:0.15em;margin-bottom:4px">Tier 1 Governing Equation — Pathway Score</div>
    Pₖ = Σⱼ (wₖⱼ · bₖⱼ) &nbsp;|&nbsp; Σⱼ wₖⱼ = 1 &nbsp;|&nbsp; Pₖ ∈ [0, 1] &nbsp;|&nbsp; Pathway activated if Pₖ ≥ {ACTIVATION_THRESHOLD}
  </div>

  <!-- Arrow -->
  <div style="display:flex;justify-content:center;align-items:center;height:32px;margin:4px 0">
    <div style="width:0;height:0;border-left:8px solid transparent;border-right:8px solid transparent;border-top:10px solid #0099cc;opacity:0.6"></div>
  </div>

  <!-- LAYER 2 -->
  <div style="font-family:'DM Mono',monospace;font-size:9px;letter-spacing:0.18em;text-transform:uppercase;color:#5a8099;text-align:center;padding:6px 0 10px">
    Layer 2 · Biological Process Estimation — Qᵢ = Σₖ∈ᵢ (vᵢₖ · Pₖ) &nbsp;|&nbsp; Σₖ vᵢₖ = 1 &nbsp;|&nbsp; Qᵢ ∈ [0, 1]
  </div>

  <div style="display:grid;grid-template-columns:1fr 1fr 1fr;gap:16px">
    {pc_inflam}
    {pc_injury}
    {pc_fibro}
  </div>

  <!-- Equation strip 2 -->
  <div style="width:100%;background:rgba(0,201,167,0.04);border:1px solid rgba(0,201,167,0.15);border-radius:10px;padding:10px 18px;margin:6px 0;font-family:'DM Mono',monospace;font-size:11px;color:#00c9a7;text-align:center">
    <div style="font-size:9px;color:#5a8099;text-transform:uppercase;letter-spacing:0.15em;margin-bottom:4px">Tier 2 & 3 Governing Equations — Aggregate + Context Adjustment</div>
    S = {PROCESS_WEIGHTS['Inflammation']}·Q_inflam + {PROCESS_WEIGHTS['Tissue Injury']}·Q_injury + {PROCESS_WEIGHTS['Fibrosis']}·Q_fibrosis = {round(S,3)} &nbsp;|&nbsp;
    S_adj = clip(S + ε, 0, 1) = {round(S_adj,3)} &nbsp;|&nbsp; RRS = S_adj × 100 = {RRS}
  </div>

  <!-- Arrow -->
  <div style="display:flex;justify-content:center;align-items:center;height:32px;margin:4px 0">
    <div style="width:0;height:0;border-left:8px solid transparent;border-right:8px solid transparent;border-top:10px solid #0099cc;opacity:0.6"></div>
  </div>

  <!-- LAYER 3 -->
  <div style="font-family:'DM Mono',monospace;font-size:9px;letter-spacing:0.18em;text-transform:uppercase;color:#5a8099;text-align:center;padding:6px 0 10px">
    Layer 3 · Clinical Output — P(rejection) = 1 / (1 + e^(−α(S_adj − θ))) = {p_rejection}
  </div>

  <div style="display:grid;grid-template-columns:1fr 1fr 1fr;gap:16px">

    <!-- RRS -->
    <div style="background:#0a1a24;border:1px solid #1a3040;border-radius:14px;padding:18px;position:relative;overflow:hidden">
      <div style="position:absolute;bottom:0;left:0;right:0;height:2px;background:linear-gradient(to right,#ff4f6e,#ff8c42)"></div>
      <div style="font-family:'DM Mono',monospace;font-size:10px;text-transform:uppercase;letter-spacing:0.12em;color:#5a8099;margin-bottom:10px">Rejection Risk Score</div>
      <div style="font-size:52px;font-weight:800;line-height:1;background:linear-gradient(135deg,#ff8c42,#ff4f6e);-webkit-background-clip:text;-webkit-text-fill-color:transparent;background-clip:text;margin-bottom:4px">{RRS}</div>
      <div style="font-family:'DM Mono',monospace;font-size:10px;color:{risk_color};text-transform:uppercase;letter-spacing:0.15em">{risk_label}</div>
      <div style="font-family:'DM Mono',monospace;font-size:9px;color:#5a8099;margin-top:6px;margin-bottom:10px">RRS = S_adj × 100 &nbsp;|&nbsp; pathway-anchored</div>
      {tier_html}
      <div style="margin-top:10px;font-family:'DM Mono',monospace;font-size:9px;color:#00c9a7;background:rgba(0,201,167,0.06);border:1px solid rgba(0,201,167,0.15);border-radius:6px;padding:5px 8px">
        P(rejection) = 1/(1+e^(−α({round(S_adj,3)}−θ))) = {p_rejection}
      </div>
    </div>

    <!-- Profile -->
    <div style="background:#0a1a24;border:1px solid #1a3040;border-radius:14px;padding:18px;position:relative;overflow:hidden">
      <div style="position:absolute;bottom:0;left:0;right:0;height:2px;background:linear-gradient(to right,#45aaf2,#a29bfe)"></div>
      <div style="font-family:'DM Mono',monospace;font-size:10px;text-transform:uppercase;letter-spacing:0.12em;color:#5a8099;margin-bottom:12px">Rejection Process Profile</div>
      {profile_html}
    </div>

    <!-- Actions -->
    <div style="background:#0a1a24;border:1px solid #1a3040;border-radius:14px;padding:18px;position:relative;overflow:hidden">
      <div style="position:absolute;bottom:0;left:0;right:0;height:2px;background:linear-gradient(to right,#00c9a7,#0099cc)"></div>
      <div style="font-family:'DM Mono',monospace;font-size:10px;text-transform:uppercase;letter-spacing:0.12em;color:#5a8099;margin-bottom:12px">Suggested Clinical Actions</div>
      {action_items(actions)}
    </div>

  </div>

  <!-- Legend -->
  <div style="display:flex;gap:20px;justify-content:center;flex-wrap:wrap;margin-top:32px;padding-top:20px;border-top:1px solid #1a3040">
    {''.join(f'<div style="display:flex;align-items:center;gap:6px;font-family:\'DM Mono\',monospace;font-size:10px;color:#5a8099"><div style="width:8px;height:8px;border-radius:50%;background:{c}"></div>{n}</div>' for n,c in [("Innate Inflammation","#ff6b6b"),("Complement","#f7b731"),("T-Cell Activation","#45aaf2"),("Antibody / Fc","#a29bfe"),("Fibrosis Remodeling","#fd9644")])}
  </div>

</div>
</body>
</html>"""


# ── STEP 6: WRITE HTML FILE AND OPEN IN BROWSER ───────────────────────────────

output_filename = "IRSM_dashboard_output.html"
with open(output_filename, "w", encoding="utf-8") as f:
    f.write(HTML)

print(f"\nDashboard saved to: {os.path.abspath(output_filename)}")
print("Opening in browser...")
webbrowser.open(f"file://{os.path.abspath(output_filename)}")
