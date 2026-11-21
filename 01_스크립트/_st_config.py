# _st_config.py - Normalises a sets_<TAG>.json model-settings file (variables, transforms, labels, rules) for 80_build_disease_panel_v.py.
#  Variable labels are read from _var_labels.json.
import os, re, json, unicodedata
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
VER = os.environ.get("ST_VER", "").strip(); SUF = f"_{VER}" if VER else ""
NFC = lambda s: unicodedata.normalize("NFC", str(s)).strip()
DISEASE_EN = {}
DISEASE_KO = {}
def _norm_transform(t):
    t = NFC(t or "none").lower()
    if t.startswith(("share_of_pop", "pct_of_pop", "ratio_pop")) or "도시지역인구" in t or "도시인구비율" in t: return "pct_of_pop"
    if t.startswith("share_of_area") or t.startswith("pct_of_area"): return "pct_of_area"
    if t in ("per100k", "per_100k", "per 100k"): return "per100k"
    if t in ("per1k", "per_1k", "per1000"): return "per1k"
    if t in ("per100k_log1p", "per100k→log1p", "per100k_log"): return "per100k_log1p"
    if t in ("log", "log1p", "none"): return t
    if "그대로" in t or t.startswith("none") or t.startswith("z"): return "none"
    if ("10만" in t or "100k" in t or "100,000" in t) and "log1p" in t: return "per100k_log1p"
    if "10만" in t or "100k" in t or "100,000" in t: return "per100k"
    if "log1p" in t: return "log1p"
    raise ValueError(f"unknown transform: {t}")
def _windows(w):
    default = {"W1_2011-15": [2011, 2015], "W2_2016-19": [2016, 2019], "W3_2020-24": [2020, 2024]}
    if not w: return default
    if isinstance(w, dict): return {k: [int(v[0]), int(v[1])] for k, v in w.items()}
    out = {}
    for x in w:
        m = re.fullmatch(r"\s*(\d{4})\s*[-~]\s*(\d{4})\s*", str(x))
        if m: out[f"W{len(out)+1}_{m.group(1)}-{m.group(2)[2:]}"] = [int(m.group(1)), int(m.group(2))]
    return out or default
def load_config(TAG):
    OUT = os.path.join(ROOT, "04_분석", TAG)
    S = json.load(open(os.path.join(OUT, f"sets_{TAG}{SUF}.json"), encoding="utf-8"))
    L = pd.DataFrame(json.load(open(os.path.join(HERE, "_var_labels.json"), encoding="utf-8"))); L["column"] = L["column"].map(NFC); L = L.set_index("column")
    hp = os.path.join(OUT, f"hypothesis_table_{TAG}.csv")
    H = pd.read_csv(hp, encoding="utf-8-sig") if os.path.exists(hp) else pd.DataFrame()
    H.columns = [c.strip() for c in H.columns]
    keycol = next((c for c in ("variable_key", "column", "variable") if c in H.columns), None)
    if keycol: H[keycol] = H[keycol].map(NFC); H = H.set_index(keycol)
    dircol = next((c for c in ("expected_direction", "예상방향") if c in H.columns), None)
    forced = [NFC(v) for v in S["forced"]]; neg = [NFC(v) for v in S.get("negative_controls", [])]
    sets = {k: [NFC(v) for v in vs] for k, vs in S["sets"].items() if not k.upper().startswith("REF_")}
    OBES = "kchs__DT_H_OBE_OBE__SR"; sens_sets = []; sens_only = []
    if TAG != "METS" and "S1_narrative_full" in sets and OBES not in sets["S1_narrative_full"]:
        sets["S1b_obesity_sensitivity"] = sets["S1_narrative_full"] + [OBES]; sens_sets.append("S1b_obesity_sensitivity")
        if not any(OBES in vs for k, vs in sets.items() if k not in sens_sets): sens_only.append(OBES)
    for _v in (S.get("sensitivity_only_extra") or []):
        if NFC(_v) not in sens_only: sens_only.append(NFC(_v))
    cand = []
    for sn, vs in sets.items():
        for v in vs:
            if v not in cand and v not in forced and v not in sens_only: cand.append(v)
    allcov = cand + [v for v in sens_only if v not in cand]
    TR = {}
    for fld in ("transform", "forced_transform", "include_transform"): TR.update({NFC(k): v for k, v in (S.get(fld) or {}).items()})
    TNOTE = str(S.get("transform_note", ""))
    def _note_rule(col):
        if not TNOTE: return "none"
        if col.startswith(("fac__", "vil__")): return "per100k_log1p" if "log1p" in TNOTE else "per100k"
        if "도시지역인구" in col: return "pct_of_pop"
        if col.startswith("chk__요양기관") or (col.startswith("nhis__") and col.endswith("진료실인원수")): return "per100k"
        return "none"
    EDIR = {NFC(k): v for k, v in (S.get("expected_direction") or {}).items()} if isinstance(S.get("expected_direction"), dict) else {}
    vars_ = {}; used = set()
    for raw in forced + allcov + neg:
        src = NFC(H.loc[raw, "source_column"]) if (keycol and raw in H.index and "source_column" in H.columns and pd.notna(H.loc[raw, "source_column"])) else raw
        tr = TR.get(raw) or (H.loc[raw, "transform"] if (keycol and raw in H.index and "transform" in H.columns and pd.notna(H.loc[raw, "transform"])) else None) or _note_rule(raw)
        tr = _norm_transform(tr if pd.notna(tr) else "none")
        if src in L.index: key, le, lk = L.loc[src, "key"], L.loc[src, "label_en"], L.loc[src, "label_ko"]
        else:
            key = re.sub(r"[^A-Za-z0-9_]+", "_", raw).strip("_").lower() or "var"
            if not re.fullmatch(r"[a-z][a-z0-9_]*", key): key = "var_" + key
            le = (H.loc[raw, "label_en"] if keycol and raw in H.index and "label_en" in H.columns else raw); lk = (H.loc[raw, "label_ko"] if keycol and raw in H.index and "label_ko" in H.columns else raw)
            print(f"[config] 라벨표에 없는 변수 → 자동 키 {key}: {raw}")
        if raw != src and raw not in L.index and "label_en" in H.columns: le, lk = H.loc[raw, "label_en"], H.loc[raw, "label_ko"]
        key = str(key); assert key not in used, f"내부 키 중복: {key} ({raw})"; used.add(key)
        ed = EDIR.get(raw) or (H.loc[raw, dircol] if (keycol and dircol and raw in H.index and pd.notna(H.loc[raw, dircol])) else "")
        vars_[raw] = dict(key=key, column=src, transform=tr, label_en=str(le), label_ko=str(lk), expected_direction=str(ed), role="forced" if raw in forced else ("negative" if raw in neg else ("candidate_sensitivity" if raw in sens_only else "candidate")))
    K = lambda raw: vars_[raw]["key"]
    RU = S.get("rules", {}) or {}
    main_set = RU.get("main_set")
    if not main_set:
        note = " ".join(str(x) for x in (RU.get("main_model_note", ""), RU.get("main_set_note", ""), S.get("main_set_policy", "")))
        hit = [sn for sn in sets if sn in note]
        main_set = hit[0] if hit and ("주모형" in note or "principal" in note.lower()) else None
    vif_action = RU.get("vif_action") or ("report" if ("유지" in " ".join(str(x) for x in (RU.get("main_model_note", ""), RU.get("main_set_note", ""))) or main_set == "S1_narrative_full") else "drop")
    rules = dict(waic_tol=float(RU.get("waic_tol", 4)), cpo_fail_frac=float(RU.get("cpo_fail_frac", 0.05)), main_set=main_set, vif_threshold=float(RU.get("vif_drop", RU.get("vif_threshold", 10))), vif_action=vif_action,
                 multiverse_n=int(RU.get("multiverse_n", RU.get("multiverse_specs", 200))), windows=_windows(RU.get("windows")),
                 covid_post_from=int(RU.get("covid_post_from", 2020)), years=RU.get("years", [2011, 2024]), raw_rules=RU)
    oc = NFC(S["outcome_column"]); m = re.match(r"^[a-z]+__(.+?)__", oc)
    off = S.get("offset") or (S.get("denominator") or {}).get("primary") or S.get("population_column") or "population"
    off = NFC(re.split(r"\s*\(", str(off))[0]); off_sens = S.get("offset_sensitivity")
    if off_sens: off_sens = NFC(re.split(r"\s*\(", str(off_sens))[0])
    off_fallback = "population_by_sex_ratio" if ("성비" in str((S.get("denominator") or {}).get("fallback", ""))) else None
    excl = list(S.get("exclude_years") or []); mm = re.search(r"(\d{4})년은[^,;]*결측", str(S.get("outcome_note", "")))
    if mm and int(mm.group(1)) not in excl: excl.append(int(mm.group(1)))
    cfg = dict(tag=TAG, outcome_column=oc, outcome_label_en=S.get("outcome_label_en", f"{DISEASE_EN.get(TAG, TAG)} treated patients (NHIS, annual)"),
               disease_ko=S.get("disease_ko", DISEASE_KO.get(TAG, m.group(1) if m else TAG)), disease_en=S.get("disease_en", DISEASE_EN.get(TAG, TAG)),
               forced=[K(v) for v in forced], candidates=[K(v) for v in cand], negative_controls=[K(v) for v in neg], sets={k: [K(v) for v in vs] for k, vs in sets.items()},
               sensitivity_sets=sens_sets, sensitivity_only_vars=[K(v) for v in sens_only],
               offset_column=off, offset_sensitivity=off_sens, offset_fallback=off_fallback, exclude_years=excl,
               vars={v["key"]: dict(v, raw_key=raw) for raw, v in vars_.items()}, rules=rules, note=S.get("note", "") or S.get("outcome_note", ""), fixed_on=S.get("fixed_on", "") or S.get("fixed_before_results", ""), source_json=f"sets_{TAG}{SUF}.json", version=VER or "v1", suffix=SUF,
               offset_E_column=S.get("offset_E_column"), offset_E_file=S.get("offset_E_file"), offset_E_sensitivity_column=S.get("offset_E_sensitivity_column"), offset_E_label=S.get("offset_E_label"))
    return cfg
if __name__ == "__main__":
    import sys; print(json.dumps(load_config(sys.argv[1]), ensure_ascii=False, indent=1)[:3000])
