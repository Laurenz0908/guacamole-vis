"""
GuaCAMOLE Diagnostics Viewer
Visualize iteration-level metrics from GuaCAMOLE's --verbose output.
"""

import json
import math
import re
from io import StringIO

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="GuaCAMOLE Diagnostics",
    page_icon="🥑",
    layout="wide",
)

# ---------------------------------------------------------------------------
# Custom CSS
# ---------------------------------------------------------------------------
st.markdown("""
<style>
    /* Tighten top margin */
    .block-container { padding-top: 1.5rem; }

    /* Metric cards */
    div[data-testid="stMetric"] {
        background: #f8f9fa;
        border: 1px solid #e9ecef;
        border-radius: 8px;
        padding: 12px 16px;
    }
    div[data-testid="stMetric"] label,
    div[data-testid="stMetric"] [data-testid="stMetricLabel"] p,
    div[data-testid="stMetricLabel"] {
        color: #495057 !important;
    }
    div[data-testid="stMetric"] [data-testid="stMetricValue"],
    div[data-testid="stMetricValue"] {
        color: #212529 !important;
    }

    /* Slider label */
    .stSlider label { font-weight: 600; }
</style>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# Helper: load iteration JSONs
# ---------------------------------------------------------------------------

def parse_iteration_files(uploaded_files):
    """Parse uploaded JSON files into a sorted list of iteration dicts.

    Handles Python-style NaN in JSON by replacing it before parsing.
    Also loads the manifest if present.
    """
    iterations = []
    manifest = None

    for f in uploaded_files:
        raw = f.read().decode("utf-8")
        name = f.name

        if name == "iteration_manifest.json":
            manifest = json.loads(raw)
            continue

        if not name.startswith("iteration_") or not name.endswith(".json"):
            continue

        # Replace NaN (not valid JSON) with null
        cleaned = raw.replace("NaN", "null")
        data = json.loads(cleaned)
        data["_filename"] = name
        iterations.append(data)

    # Sort by (iteration, rep)
    iterations.sort(key=lambda d: (d["iteration"], d["rep"]))
    return iterations, manifest


def build_step_labels(iterations):
    """Build human-readable labels for each iteration step."""
    labels = []
    for it in iterations:
        label = f"Cycle {it['iteration']}  /  Rep {it['rep']}  (T={it['threshold']:.3g})"
        labels.append(label)
    return labels


def taxid_name_map(guac_file):
    """Parse a .guac output file to get taxid -> species name mapping."""
    mapping = {}
    content = guac_file.read().decode("utf-8")
    for line in content.strip().split("\n")[1:]:  # skip header
        parts = line.split("\t")
        if len(parts) >= 2:
            name, taxid = parts[0], parts[1]
            mapping[int(taxid)] = name
    return mapping


# ---------------------------------------------------------------------------
# Sidebar: File upload
# ---------------------------------------------------------------------------
st.sidebar.title("🥑 GuaCAMOLE Diagnostics")
st.sidebar.markdown("Upload the `--verbose` output files from a GuaCAMOLE run.")

uploaded_jsons = st.sidebar.file_uploader(
    "Iteration JSON files + manifest",
    type=["json"],
    accept_multiple_files=True,
    help="Select all `iteration_*.json` files and `iteration_manifest.json`",
)

uploaded_guac = st.sidebar.file_uploader(
    "GuaCAMOLE output file (.guac) — optional",
    type=["guac", "txt", "tsv"],
    help="Provides species names for taxids",
)

# ---------------------------------------------------------------------------
# Main: require data
# ---------------------------------------------------------------------------
if not uploaded_jsons:
    st.title("🥑 GuaCAMOLE Diagnostics Viewer")
    st.info(
        "Upload iteration JSON files from a GuaCAMOLE `--verbose` run using the sidebar.\n\n"
        "Optionally upload the `.guac` output file to display species names instead of taxids."
    )
    st.stop()

# ---------------------------------------------------------------------------
# Parse data
# ---------------------------------------------------------------------------
iterations, manifest = parse_iteration_files(uploaded_jsons)

if len(iterations) == 0:
    st.error("No valid iteration JSON files found. Files must be named `iteration_<N>_rep_<M>.json`.")
    st.stop()

# Taxid -> name mapping
name_map = {}
if uploaded_guac:
    name_map = taxid_name_map(uploaded_guac)


def taxid_label(tid):
    """Return a readable label for a taxid."""
    if tid in name_map:
        return f"{name_map[tid]} ({tid})"
    return str(tid)


# ---------------------------------------------------------------------------
# Iteration slider
# ---------------------------------------------------------------------------
st.title("🥑 GuaCAMOLE Diagnostics Viewer")

step_labels = build_step_labels(iterations)
n_steps = len(iterations)

if n_steps == 1:
    step_idx = 0
    st.markdown(f"**Step:** {step_labels[0]}")
else:
    step_idx = st.slider(
        "Iteration step",
        min_value=0,
        max_value=n_steps - 1,
        value=0,
        format=f"Step %d of {n_steps}",
    )
    st.markdown(f"**{step_labels[step_idx]}**")

current = iterations[step_idx]

# ---------------------------------------------------------------------------
# Summary metrics
# ---------------------------------------------------------------------------
col1, col2, col3, col4 = st.columns(4)
col1.metric("Cycle", current["iteration"])
col2.metric("Repetition", current["rep"])
col3.metric("Threshold", f"{current['threshold']:.4g}")
col4.metric("Taxa remaining", len(current["taxids"]))

n_removed = len(current.get("removed_taxids", []))
if n_removed > 0:
    removed_names = [taxid_label(t) for t in current["removed_taxids"]]
    st.warning(f"**{n_removed} taxa removed this step:** {', '.join(removed_names)}")

# ---------------------------------------------------------------------------
# Tabs for the different visualizations
# ---------------------------------------------------------------------------
tab_eff, tab_res, tab_ab, tab_timeline = st.tabs([
    "📈 Sequencing Efficiency",
    "📊 Residuals",
    "📉 Abundances",
    "🗓️ Taxa Removal Timeline",
])

# ===== TAB 1: Efficiency curve =============================================
with tab_eff:
    eff = np.array(current["efficiencies"])
    gc_bins = np.arange(len(eff))

    # Mask zero-regions (no data at extreme GC)
    eff_masked = eff.copy()
    eff_masked[eff_masked == 0] = np.nan

    fig_eff = go.Figure()
    fig_eff.add_trace(go.Scatter(
        x=gc_bins,
        y=eff_masked,
        mode="lines",
        line=dict(color="#2a9d8f", width=2.5),
        name="Efficiency",
        hovertemplate="GC bin: %{x}<br>Efficiency: %{y:.4f}<extra></extra>",
    ))

    fig_eff.update_layout(
        title="GC-dependent Sequencing Efficiency",
        xaxis_title="GC bin (%)",
        yaxis_title="Relative Sequencing Efficiency",
        yaxis_range=[0, 1.1],
        template="plotly_white",
        height=450,
        margin=dict(t=50, b=50),
    )

    # Show comparison with previous step
    if step_idx > 0:
        prev_eff = np.array(iterations[step_idx - 1]["efficiencies"])
        prev_eff[prev_eff == 0] = np.nan
        fig_eff.add_trace(go.Scatter(
            x=gc_bins,
            y=prev_eff,
            mode="lines",
            line=dict(color="#adb5bd", width=1.5, dash="dot"),
            name="Previous step",
            hovertemplate="GC bin: %{x}<br>Efficiency: %{y:.4f}<extra></extra>",
        ))

    st.plotly_chart(fig_eff, use_container_width=True)

# ===== TAB 2: Residuals ====================================================
with tab_res:
    residuals = np.array(current["residuals"])  # shape: (n_gc_bins, n_taxa)
    taxids = current["taxids"]
    threshold = current["threshold"]
    res_range = np.array(current["res_range"])

    # Identify taxa that will be removed (res_range > threshold)
    will_remove = set()
    for i, rr in enumerate(res_range):
        if rr > threshold:
            will_remove.add(i)

    # Taxa selection
    all_taxid_labels = [taxid_label(t) for t in taxids]

    # Default: show taxa to be removed + a few top-residual taxa
    default_indices = sorted(will_remove)
    if len(default_indices) < 5:
        # Add some high res_range taxa
        sorted_by_res = np.argsort(res_range)[::-1]
        for idx in sorted_by_res:
            if idx not in will_remove:
                default_indices.append(idx)
            if len(default_indices) >= 8:
                break

    default_selection = [all_taxid_labels[i] for i in default_indices if i < len(all_taxid_labels)]

    selected_taxa = st.multiselect(
        "Select taxa to display",
        options=all_taxid_labels,
        default=default_selection,
        help="Red lines = taxa exceeding residual range threshold (will be removed)",
    )

    selected_indices = [all_taxid_labels.index(s) for s in selected_taxa]

    gc_bins_res = np.arange(residuals.shape[0])

    fig_res = go.Figure()

    for idx in selected_indices:
        tid = taxids[idx]
        label = taxid_label(tid)
        resid_col = residuals[:, idx]
        is_flagged = idx in will_remove

        fig_res.add_trace(go.Scatter(
            x=gc_bins_res,
            y=resid_col,
            mode="lines",
            name=f"{'🔴 ' if is_flagged else ''}{label}",
            line=dict(
                color="crimson" if is_flagged else None,
                width=2.5 if is_flagged else 1.5,
            ),
            opacity=1.0 if is_flagged else 0.7,
            hovertemplate=(
                f"{label}<br>"
                "GC bin: %{x}<br>"
                "Residual: %{y:.4f}<br>"
                f"Range: {res_range[idx]:.4f}"
                "<extra></extra>"
            ),
        ))

    # Add threshold band
    fig_res.add_hline(
        y=threshold / 2, line_dash="dash", line_color="rgba(220,53,69,0.3)",
        annotation_text=f"±T/2 ({threshold/2:.3g})",
        annotation_position="top right",
    )
    fig_res.add_hline(
        y=-threshold / 2, line_dash="dash", line_color="rgba(220,53,69,0.3)",
    )

    fig_res.update_layout(
        title="Relative Residuals per Taxon",
        xaxis_title="GC bin (%)",
        yaxis_title="Relative Residual (Obs − Exp) / Exp",
        template="plotly_white",
        height=550,
        margin=dict(t=50, b=50),
        legend=dict(
            font=dict(size=10),
            yanchor="top",
            y=-0.15,
            xanchor="left",
            x=0,
            orientation="h",
        ),
    )

    st.plotly_chart(fig_res, use_container_width=True)

    # Residual range bar chart
    st.markdown("#### Residual Range per Taxon")

    res_df = pd.DataFrame({
        "taxid": taxids,
        "label": all_taxid_labels,
        "res_range": res_range,
        "flagged": [i in will_remove for i in range(len(taxids))],
    }).sort_values("res_range", ascending=True)

    colors = ["crimson" if f else "#2a9d8f" for f in res_df["flagged"]]

    fig_bar = go.Figure()
    fig_bar.add_trace(go.Bar(
        y=res_df["label"],
        x=res_df["res_range"],
        orientation="h",
        marker_color=colors,
        hovertemplate="%{y}<br>Range: %{x:.4f}<extra></extra>",
    ))

    fig_bar.add_vline(
        x=threshold, line_dash="dash", line_color="crimson",
        annotation_text=f"Threshold ({threshold:.3g})",
    )

    fig_bar.update_layout(
        template="plotly_white",
        height=max(400, len(taxids) * 18),
        margin=dict(l=200, t=30, b=40),
        xaxis_title="Residual Range (max − min)",
    )

    st.plotly_chart(fig_bar, use_container_width=True)


# ===== TAB 3: Abundances ===================================================
with tab_ab:
    ab = np.array(current["abundances"])
    taxids_ab = current["taxids"]
    labels_ab = [taxid_label(t) for t in taxids_ab]

    ab_df = pd.DataFrame({
        "label": labels_ab,
        "abundance": ab,
    }).sort_values("abundance", ascending=True)

    fig_ab = go.Figure()
    fig_ab.add_trace(go.Bar(
        y=ab_df["label"],
        x=ab_df["abundance"],
        orientation="h",
        marker_color="#264653",
        hovertemplate="%{y}<br>Abundance: %{x:.6f}<extra></extra>",
    ))

    fig_ab.update_layout(
        title=f"Abundance Estimates — Cycle {current['iteration']}, Rep {current['rep']}",
        xaxis_title="Relative Abundance",
        template="plotly_white",
        height=max(400, len(taxids_ab) * 18),
        margin=dict(l=200, t=50, b=40),
    )

    st.plotly_chart(fig_ab, use_container_width=True)

    # Abundance change compared to previous step
    if step_idx > 0:
        st.markdown("#### Abundance Change vs. Previous Step")
        prev = iterations[step_idx - 1]
        prev_ab_dict = dict(zip(prev["taxids"], prev["abundances"]))

        changes = []
        for tid, ab_val in zip(taxids_ab, ab):
            if tid in prev_ab_dict and prev_ab_dict[tid] > 0:
                fold_change = ab_val / prev_ab_dict[tid]
                changes.append({
                    "label": taxid_label(tid),
                    "log2_fc": np.log2(fold_change) if fold_change > 0 else 0,
                })

        if changes:
            change_df = pd.DataFrame(changes).sort_values("log2_fc")
            colors_fc = ["#e76f51" if v > 0 else "#264653" for v in change_df["log2_fc"]]

            fig_fc = go.Figure()
            fig_fc.add_trace(go.Bar(
                y=change_df["label"],
                x=change_df["log2_fc"],
                orientation="h",
                marker_color=colors_fc,
                hovertemplate="%{y}<br>log₂ FC: %{x:.4f}<extra></extra>",
            ))
            fig_fc.update_layout(
                xaxis_title="log₂ Fold Change",
                template="plotly_white",
                height=max(400, len(changes) * 18),
                margin=dict(l=200, t=30, b=40),
            )
            fig_fc.add_vline(x=0, line_color="gray", line_width=1)
            st.plotly_chart(fig_fc, use_container_width=True)


# ===== TAB 4: Taxa Removal Timeline ========================================
with tab_timeline:
    # Build a table: for each iteration step, which taxa were removed
    removal_events = []
    cumulative_removed = set()

    for i, it in enumerate(iterations):
        removed = it.get("removed_taxids", [])
        for tid in removed:
            removal_events.append({
                "Taxon": taxid_label(tid),
                "Taxid": tid,
                "Removed at Cycle": it["iteration"],
                "Removed at Rep": it["rep"],
                "Threshold": it["threshold"],
            })
            cumulative_removed.add(tid)

    if manifest:
        # Check for taxa in manifest but not yet seen (removed before first iteration?)
        for tid in manifest.get("all_skipped_taxids", []):
            if tid not in cumulative_removed:
                removal_events.append({
                    "Taxon": taxid_label(tid),
                    "Taxid": tid,
                    "Removed at Cycle": "?",
                    "Removed at Rep": "?",
                    "Threshold": "—",
                })

    if not removal_events:
        st.info("No taxa were removed during the iterations.")
    else:
        removal_df = pd.DataFrame(removal_events)

        # Summary chart: number of taxa removed per cycle
        if all(isinstance(c, (int, float)) for c in removal_df["Removed at Cycle"]):
            cycle_counts = removal_df.groupby("Removed at Cycle").size().reset_index(name="Count")

            fig_timeline = go.Figure()
            fig_timeline.add_trace(go.Bar(
                x=cycle_counts["Removed at Cycle"],
                y=cycle_counts["Count"],
                marker_color="#e76f51",
                text=cycle_counts["Count"],
                textposition="outside",
            ))
            fig_timeline.update_layout(
                title="Taxa Removed per Cycle",
                xaxis_title="Cycle",
                yaxis_title="Number of taxa removed",
                template="plotly_white",
                height=350,
                margin=dict(t=50, b=40),
                xaxis=dict(dtick=1),
            )
            st.plotly_chart(fig_timeline, use_container_width=True)

        # Taxa remaining over time
        st.markdown("#### Taxa Count Over Iterations")
        taxa_counts = [{"step": build_step_labels([it])[0], "n_taxa": len(it["taxids"])} for it in iterations]
        taxa_count_df = pd.DataFrame(taxa_counts)

        fig_count = go.Figure()
        fig_count.add_trace(go.Scatter(
            x=list(range(len(taxa_count_df))),
            y=taxa_count_df["n_taxa"],
            mode="lines+markers",
            line=dict(color="#2a9d8f", width=2.5),
            marker=dict(size=8),
            hovertemplate="%{text}<br>Taxa: %{y}<extra></extra>",
            text=taxa_count_df["step"],
        ))
        fig_count.update_layout(
            xaxis_title="Step",
            yaxis_title="Number of taxa",
            template="plotly_white",
            height=300,
            margin=dict(t=30, b=40),
        )
        st.plotly_chart(fig_count, use_container_width=True)

        # Detailed table
        st.markdown("#### Removal Details")
        st.dataframe(removal_df, use_container_width=True, hide_index=True)