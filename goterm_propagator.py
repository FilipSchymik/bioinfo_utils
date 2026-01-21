from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd


GO_ROOTS = {"bp": "GO:0008150", "mf": "GO:0003674", "cc": "GO:0005575"}
ROOT_IDS = set(GO_ROOTS.values())


@dataclass(frozen=True)
class GOTerm:
    go_id: str
    name: str
    namespace: str
    parents_is_a: Tuple[str, ...]
    parents_part_of: Tuple[str, ...]
    is_obsolete: bool
    replaced_by: Tuple[str, ...]
    consider: Tuple[str, ...]


def _ns_to_aspect(ns: str) -> Optional[str]:
    return {
        "biological_process": "bp",
        "molecular_function": "mf",
        "cellular_component": "cc",
    }.get(ns)


def parse_go_obo(obo_path: str) -> Dict[str, GOTerm]:
    terms: Dict[str, GOTerm] = {}
    current: Dict[str, object] = {}
    in_term = False

    def flush():
        nonlocal current, in_term
        if not in_term:
            return
        go_id = current.get("id")
        if not go_id:
            current = {}
            in_term = False
            return

        terms[str(go_id)] = GOTerm(
            go_id=str(go_id),
            name=str(current.get("name", "")),
            namespace=str(current.get("namespace", "")),
            parents_is_a=tuple(current.get("is_a", ())),
            parents_part_of=tuple(current.get("part_of", ())),
            is_obsolete=bool(current.get("is_obsolete", False)),
            replaced_by=tuple(current.get("replaced_by", ())),
            consider=tuple(current.get("consider", ())),
        )
        current = {}
        in_term = False

    with open(obo_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line == "[Term]":
                flush()
                in_term = True
                current = {"is_a": [], "part_of": [], "replaced_by": [], "consider": []}
                continue
            if not in_term:
                continue

            if line.startswith("id: "):
                current["id"] = line.split("id: ", 1)[1].strip()
            elif line.startswith("name: "):
                current["name"] = line.split("name: ", 1)[1].strip()
            elif line.startswith("namespace: "):
                current["namespace"] = line.split("namespace: ", 1)[1].strip()
            elif line.startswith("is_obsolete: "):
                current["is_obsolete"] = (line.split("is_obsolete: ", 1)[1].strip().lower() == "true")
            elif line.startswith("is_a: "):
                parent = line.split("is_a: ", 1)[1].split("!", 1)[0].strip()
                current["is_a"].append(parent)
            elif line.startswith("relationship: "):
                rest = line.split("relationship: ", 1)[1].strip()
                if rest.startswith("part_of "):
                    parent = rest.split("part_of ", 1)[1].split("!", 1)[0].strip()
                    current["part_of"].append(parent)
            elif line.startswith("replaced_by: "):
                rid = line.split("replaced_by: ", 1)[1].strip()
                current["replaced_by"].append(rid)
            elif line.startswith("consider: "):
                cid = line.split("consider: ", 1)[1].strip()
                current["consider"].append(cid)

    flush()
    return terms


def propagate_go_annotations(
    df: pd.DataFrame,
    obo_path: str,
    protein_col: str = "Protein",
    term_col: str = "GO_term/EC",
    score_col: str = "Score",
    ic_df: Optional[pd.DataFrame] = None,
    ic_term_col: str = "go_term",
    ic_value_col: str = "IC",
    keep_original_terms: bool = True,
    relations: Tuple[str, ...] = ("is_a",),
    exclude_roots: bool = True,
    obsolete_mode: str = "drop",           # "drop" | "map"
    obsolete_prefer: str = "replaced_by",  # "replaced_by" | "consider"
    score_min: Optional[float] = None,     # NEW: filter before propagation
    round_decimals: int = 3,               # NEW: rounding
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns (propagated_df, stats_df).

    - Filters input annotations by score_min (if provided) BEFORE propagation.
    - Rounds Score/IC in propagated_df and percent in stats_df to round_decimals.

    IC is joined only by exact GO ids from ic_df (never inherited).
    """
    go = parse_go_obo(obo_path)

    use_part_of = "part_of" in relations
    use_is_a = "is_a" in relations
    if not use_is_a and not use_part_of:
        raise ValueError("relations must include at least one of: 'is_a', 'part_of'")
    if obsolete_mode not in {"drop", "map"}:
        raise ValueError("obsolete_mode must be 'drop' or 'map'")
    if obsolete_prefer not in {"replaced_by", "consider"}:
        raise ValueError("obsolete_prefer must be 'replaced_by' or 'consider'")

    # strict IC lookup
    ic_map: Dict[str, float] = {}
    if ic_df is not None:
        if ic_term_col not in ic_df.columns or ic_value_col not in ic_df.columns:
            raise ValueError(f"ic_df must have columns: {ic_term_col!r}, {ic_value_col!r}")
        _ic = ic_df[[ic_term_col, ic_value_col]].dropna()
        _ic[ic_term_col] = _ic[ic_term_col].astype(str)
        ic_map = _ic.groupby(ic_term_col)[ic_value_col].max().to_dict()

    @lru_cache(maxsize=None)
    def ancestors(term_id: str) -> Tuple[str, ...]:
        """Get all ancestor terms for a given term, excluding roots if requested."""
        term = go.get(term_id)
        if term is None or term.is_obsolete:
            return tuple()

        parents: List[str] = []
        if use_is_a:
            parents.extend(term.parents_is_a)
        if use_part_of:
            parents.extend(term.parents_part_of)

        if not parents:
            return tuple()

        out: List[str] = []
        stack = list(parents)
        seen: Set[str] = set(parents)

        while stack:
            parent_id = stack.pop()
            
            # Skip if already processed or is root
            if exclude_roots and parent_id in ROOT_IDS:
                continue
            
            parent_term = go.get(parent_id)
            if parent_term is None or parent_term.is_obsolete:
                continue

            out.append(parent_id)

            # Add grandparents to stack
            grandparent_ids: List[str] = []
            if use_is_a:
                grandparent_ids.extend(parent_term.parents_is_a)
            if use_part_of:
                grandparent_ids.extend(parent_term.parents_part_of)
            
            for gp_id in grandparent_ids:
                if gp_id not in seen:
                    seen.add(gp_id)
                    stack.append(gp_id)

        return tuple(out)

    @lru_cache(maxsize=None)
    def _resolve_term(term_id: str) -> Tuple[str, ...]:
        """Resolve obsolete terms to their replacements if mapping is enabled."""
        term = go.get(term_id)
        if term is None:
            return tuple()
        
        if not term.is_obsolete:
            return (term_id,)

        if obsolete_mode == "drop":
            return tuple()

        # Get replacement candidates (primary preference, then secondary)
        primary = term.replaced_by if obsolete_prefer == "replaced_by" else term.consider
        secondary = term.consider if obsolete_prefer == "replaced_by" else term.replaced_by
        candidates = list(primary) + list(secondary)
        
        # Return only valid, non-obsolete replacements
        return tuple(
            c for c in candidates 
            if (c in go and not go[c].is_obsolete)
        )

    # ---- Input filtering and stats collection
    n_input_rows = int(df.shape[0])
    
    # Extract and filter required columns
    sub = df[[protein_col, term_col, score_col]].dropna(
        subset=[protein_col, term_col, score_col]
    ).copy()
    n_rows_nonnull = int(sub.shape[0])

    # Apply score filter if specified
    n_rows_dropped_score = 0
    if score_min is not None:
        sub[score_col] = pd.to_numeric(sub[score_col], errors="coerce")
        n_before = int(sub.shape[0])
        sub = sub.dropna(subset=[score_col])
        sub = sub[sub[score_col] >= float(score_min)]
        n_after = int(sub.shape[0])
        n_rows_dropped_score = n_before - n_after

    sub[term_col] = sub[term_col].astype(str)
    input_unique_terms = int(sub[term_col].nunique())

    # Filter GO-like terms and compute stats in a single pass
    is_go_like = sub[term_col].str.startswith("GO:")
    n_rows_go_like = int(is_go_like.sum())
    
    # Use vectorized operations for better performance
    go_like_mask = is_go_like
    if go_like_mask.any():
        go_like_terms_series = sub.loc[go_like_mask, term_col]
        # Vectorized check: use isin() for better performance
        known_mask = go_like_terms_series.isin(go)
        n_rows_known = int(known_mask.sum())
        n_rows_unknown = int((~known_mask).sum())
        
        # Count unique unknown GO terms
        unknown_terms_series = go_like_terms_series[~known_mask]
        n_unique_unknown_terms = int(unknown_terms_series.nunique()) if n_rows_unknown > 0 else 0
        
        # For obsolete check, use vectorized operations
        known_terms_series = go_like_terms_series[known_mask]
        if len(known_terms_series) > 0:
            # Map to is_obsolete using vectorized apply (faster than lambda map)
            obsolete_mask = known_terms_series.map(lambda tid: go[tid].is_obsolete)
            n_rows_obsolete = int(obsolete_mask.sum())
            n_rows_nonobsolete = int((~obsolete_mask).sum())
            
            # Count unique obsolete GO terms
            obsolete_terms_series = known_terms_series[obsolete_mask]
            n_unique_obsolete_terms = int(obsolete_terms_series.nunique()) if n_rows_obsolete > 0 else 0
            
            # Count mappable obsolete terms
            n_rows_obsolete_mapped = 0
            if obsolete_mode == "map" and n_rows_obsolete > 0:
                obsolete_terms = obsolete_terms_series.tolist()
                # Cache results to avoid redundant calls
                resolved_cache = {tid: _resolve_term(tid) for tid in obsolete_terms}
                n_rows_obsolete_mapped = sum(1 for tid in obsolete_terms if len(resolved_cache[tid]) > 0)
        else:
            # No known terms, so all GO-like terms are unknown
            n_rows_obsolete = 0
            n_rows_nonobsolete = 0
            n_rows_obsolete_mapped = 0
            n_unique_obsolete_terms = 0
    else:
        n_rows_known = 0
        n_rows_unknown = 0
        n_rows_obsolete = 0
        n_rows_nonobsolete = 0
        n_rows_obsolete_mapped = 0
        n_unique_unknown_terms = 0
        n_unique_obsolete_terms = 0

    # ---- Propagation accumulation
    max_score: Dict[Tuple[str, str], float] = {}
    propagated_flag: Dict[Tuple[str, str], bool] = {}

    # Pre-convert score column to float for better performance
    if score_min is None:
        sub[score_col] = pd.to_numeric(sub[score_col], errors="coerce")
    
    # Single pass through data for propagation
    for prot, raw_tid, score in sub.itertuples(index=False, name=None):
        raw_tid = str(raw_tid)
        if not raw_tid.startswith("GO:"):
            continue

        seeds = _resolve_term(raw_tid)
        if not seeds:
            continue

        score_float = float(score)
        
        # Process seed terms (original terms)
        for seed_term_id in seeds:
            if keep_original_terms and not (exclude_roots and seed_term_id in ROOT_IDS):
                key = (prot, seed_term_id)
                # Update max score if needed
                current_max = max_score.get(key)
                if current_max is None or score_float > current_max:
                    max_score[key] = score_float
                propagated_flag[key] = False

            # Propagate to ancestors
            for ancestor_id in ancestors(seed_term_id):
                key = (prot, ancestor_id)
                # Update max score if needed
                current_max = max_score.get(key)
                if current_max is None or score_float > current_max:
                    max_score[key] = score_float
                # Mark as propagated (only if not already set to False)
                if key not in propagated_flag:
                    propagated_flag[key] = True

    # ---- Build output dataframe
    out_rows = []
    # Pre-compute term lookups to avoid repeated dictionary access
    for (prot, term_id), score in max_score.items():
        term = go[term_id]
        out_rows.append(
            (
                prot,
                term_id,
                score,
                term.name,
                _ns_to_aspect(term.namespace),
                ic_map.get(term_id, float("nan")),
                propagated_flag.get((prot, term_id), True),
            )
        )

    out = pd.DataFrame(
        out_rows, 
        columns=["Protein", "GO_term", "Score", "Name", "Aspect", "IC", "Propagated"]
    )
    out.sort_values(
        ["Protein", "Propagated", "Score", "GO_term"], 
        ascending=[True, True, False, True], 
        inplace=True
    )
    out.reset_index(drop=True, inplace=True)

    # ---- Stats derived from output
    n_output_rows = int(out.shape[0])
    n_output_unique_terms = int(out["GO_term"].nunique()) if n_output_rows > 0 else 0
    n_output_original_rows = int((~out["Propagated"]).sum()) if n_output_rows > 0 else 0
    n_output_propagated_rows = int(out["Propagated"].sum()) if n_output_rows > 0 else 0
    new_unique_terms = int(out.loc[out["Propagated"], "GO_term"].nunique()) if n_output_rows > 0 else 0

    def pct(x: int, denom: int) -> float:
        return (100.0 * x / denom) if denom > 0 else 0.0

    # Build stats dataframe
    stats = [
        ("input_rows_total", n_input_rows, 100.0),
        ("input_rows_nonnull_cols", n_rows_nonnull, pct(n_rows_nonnull, n_input_rows)),
    ]
    if score_min is not None:
        stats.extend([
            ("input_rows_dropped_by_score", n_rows_dropped_score, pct(n_rows_dropped_score, n_rows_nonnull)),
            ("input_rows_after_score_filter", int(sub.shape[0]), pct(int(sub.shape[0]), n_rows_nonnull)),
        ])

    stats.extend([
        ("input_unique_terms_raw_after_filters", input_unique_terms, 100.0),
        ("input_rows_GO_like", n_rows_go_like, pct(n_rows_go_like, int(sub.shape[0]))),
        ("input_rows_GO_known", n_rows_known, pct(n_rows_known, n_rows_go_like) if n_rows_go_like > 0 else 0.0),
        ("input_rows_GO_unknown_dropped", n_rows_unknown, pct(n_rows_unknown, n_rows_go_like) if n_rows_go_like > 0 else 0.0),
        ("input_unique_GO_terms_unknown", n_unique_unknown_terms, pct(n_unique_unknown_terms, input_unique_terms) if input_unique_terms > 0 else 0.0),
        ("input_rows_GO_obsolete", n_rows_obsolete, pct(n_rows_obsolete, n_rows_known) if n_rows_known > 0 else 0.0),
        ("input_unique_GO_terms_obsolete", n_unique_obsolete_terms, pct(n_unique_obsolete_terms, input_unique_terms) if input_unique_terms > 0 else 0.0),
        ("input_rows_GO_nonobsolete_used", n_rows_nonobsolete, pct(n_rows_nonobsolete, n_rows_known) if n_rows_known > 0 else 0.0),
    ])
    
    if obsolete_mode == "map":
        stats.append(
            ("input_rows_GO_obsolete_mapped_to_replacement", n_rows_obsolete_mapped, 
             pct(n_rows_obsolete_mapped, n_rows_obsolete) if n_rows_obsolete > 0 else 0.0)
        )

    stats.extend([
        ("output_rows_total", n_output_rows, 100.0),
        ("output_unique_terms", n_output_unique_terms, 100.0 if n_output_unique_terms > 0 else 0.0),
        ("output_rows_original", n_output_original_rows, pct(n_output_original_rows, n_output_rows) if n_output_rows > 0 else 0.0),
        ("output_rows_propagated", n_output_propagated_rows, pct(n_output_propagated_rows, n_output_rows) if n_output_rows > 0 else 0.0),
        ("output_unique_terms_new_from_propagation", new_unique_terms, pct(new_unique_terms, n_output_unique_terms) if n_output_unique_terms > 0 else 0.0),
    ])

    stats_df = pd.DataFrame(stats, columns=["metric", "count", "percent"])

    # ---- Rounding
    if round_decimals is not None:
        out["Score"] = out["Score"].round(round_decimals)
        out["IC"] = out["IC"].round(round_decimals)
        stats_df["percent"] = stats_df["percent"].round(round_decimals)

    return out, stats_df
