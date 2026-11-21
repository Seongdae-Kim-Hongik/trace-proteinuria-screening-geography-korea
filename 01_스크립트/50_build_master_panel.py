# 50_build_master_panel.py - Merges the district-year source tables in 02_패널/ into one district x year master panel and a variable dictionary.
#  Input: 02_패널/*_wide.csv, district master key, population panel (ST_DATA_ROOT/식중독최종.csv).
#  Output: 02_패널/MASTER_시군구_연도_wide.parquet
import os, re, sys, importlib.util, datetime
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAN = os.path.join(ROOT, "02_패널")
FBD = os.path.join(os.environ["ST_DATA_ROOT"], "식중독최종.csv")

_spec = importlib.util.spec_from_file_location("region_crosswalk", os.path.join(ROOT, "01_스크립트", "40_region_crosswalk.py"))
rc = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(rc)
sido_std, sgg_clean, make_key = rc.sido_std, rc.sgg_clean, rc.make_key

RENAME = {
    ("경기도", "여주군"): ("경기도", "여주시", "개명"), ("충청남도", "당진군"): ("충청남도", "당진시", "개명"),
    ("경기도", "양주군"): ("경기도", "양주시", "개명"), ("경기도", "포천군"): ("경기도", "포천시", "개명"),
    ("경기도", "광주군"): ("경기도", "광주시", "개명"), ("경기도", "화성군"): ("경기도", "화성시", "개명"),
    ("경상남도", "통합창원시"): ("경상남도", "창원시", ""), ("충청북도", "통합청주시"): ("충청북도", "청주시", ""),
    ("세종시", "세종"): ("세종시", "세종시", ""), ("세종시", "세종특별자치시"): ("세종시", "세종시", ""),
    ("울산시", "울주"): ("울산시", "울주군", "표기정정"), ("인천시", "강화"): ("인천시", "강화군", "표기정정"),
    ("경기도", "의정부"): ("경기도", "의정부시", "표기정정"), ("경기도", "얀양시"): ("경기도", "안양시", "표기정정"),
    ("경기도", "수원사"): ("경기도", "수원시", "표기정정"), ("강원도", "회천군"): ("강원도", "화천군", "표기정정"),
    ("경상북도", "당진시"): ("충청남도", "당진시", "표기정정"), ("서울시", "광명시"): ("경기도", "광명시", "표기정정"),
    ("서울시", "인천계양구"): ("인천시", "계양구", "표기정정"), ("서울시", "고양시덕양구"): ("경기도", "고양시덕양구", "표기정정"),
    ("충청북도", "상당구"): ("충청북도", "청주시상당구", "표기정정"), ("충청북도", "흥덕구"): ("충청북도", "청주시흥덕구", "표기정정"),
}
MERGE = {
    ("충청북도", "청원군"): ("충청북도", "청주시", "청원군합산", "keep"),
    ("경상남도", "마산시"): ("경상남도", "창원시", "창원통합", "nan"),
    ("경상남도", "진해시"): ("경상남도", "창원시", "창원통합", "nan"),
    ("제주도", "북제주군"): ("제주도", "제주시", "북제주군합산", "keep"),
    ("제주도", "남제주군"): ("제주도", "서귀포시", "남제주군합산", "keep"),
}
PARENT_CITIES = ["고양시", "성남시", "수원시", "안산시", "안양시", "용인시", "부천시", "청주시", "천안시", "전주시", "포항시", "창원시", "화성시"]


def map_region(sido_raw, sgg_raw, master_keys, subs_unmatched=False, merged_unmatched=False):
    sgg_raw = re.sub(r"\([^)]*~[^)]*\)|\(변동전\)", "", str(sgg_raw))
    k = make_key(sido_raw, sgg_raw)
    sd, sg = k.split("|", 1)
    role, flag = "primary", ""
    if (sd, sg) in RENAME:
        sd, sg, flag = RENAME[(sd, sg)]
    ratio_mode = None
    if (sd, sg) in MERGE:
        if merged_unmatched:
            return dict(key=None, sido_std=sd, sigungu=sg, role="merged", flag="폐지지역(미매칭)", ratio_mode=None)
        sd, sg, flag, ratio_mode = MERGE[(sd, sg)]; role = "merged"
    else:
        for city in PARENT_CITIES:
            if sg.startswith(city) and sg != city and sg.endswith("구"):
                if subs_unmatched:
                    return dict(key=None, sido_std=sd, sigungu=sg, role="sub", flag="하위구(미매칭)", ratio_mode=None)
                sd, sg, role, flag = sd, city, "sub", "하위구합산"; ratio_mode = "sub"; break
    key = f"{sd}|{sg}"
    if key not in master_keys:
        return dict(key=None, sido_std=sd, sigungu=sg, role=role, flag=flag, ratio_mode=ratio_mode)
    return dict(key=key, sido_std=sd, sigungu=sg, role=role, flag=flag, ratio_mode=ratio_mode)


def _exists(f): return os.path.exists(os.path.join(PAN, f))

ID_DROP = ["period", "sgg_code", "sgg_code_std", "level", "parent_code", "code_note", "region_nm", "kchs_code",
           "sido_clean", "sigungu_clean", "sido_short", "sheet", "source_file"]

SOURCES = [
    dict(name="kchs", prefix="kchs__", files=["KCHS_시군구_연도_wide.csv"], sido="sido", sgg="sigungu", merge_counts=False,
         ratio=lambda c: not c.endswith("__N"), label="지역사회건강조사(KCHS)", pool="kchs"),
    dict(name="nhis", prefix="nhis__", files=["NHIS_주요질환_시군구_연도_wide.csv", "NHIS_진료현황_시군구_연도_wide.csv"], sido="sido", sgg="sigungu",
         merge_counts=True, ratio=lambda c: False, label="건보 지역별의료이용통계(NHIS)"),
    dict(name="chk", prefix="chk__", files=["NHIS검진_시군구_연도_wide.csv"], sido="sido", sgg="sigungu", merge_counts=True,
         ratio="dict", label="건강검진통계(NHIS 검진)"),
    dict(name="mort", prefix="mort__", files=["사망원인_시군구_연도_wide.csv"], sido="sido", sgg="sigungu", merge_counts=True, level="시군구",
         ratio=lambda c: c.endswith("_rate") or c.endswith("_asr"), label="사망원인통계"),
    dict(name="pop", prefix="pop__", files=["인구동향_시군구_연도_wide.csv"], sido="sido", sgg="sigungu", merge_counts=True, level="시군구",
         ratio=lambda c: ("출산율" in c) or ("평균" in c) or ("연령" in c and "출생아수" not in c), label="인구동향조사"),
    dict(name="fac", prefix="fac__", files=["체육시설_시군구_연도_wide.csv"], sido="sido_clean", sgg="sigungu_clean", merge_counts=True,
         ratio=lambda c: False, label="문체부 전국공공체육시설현황(공공체육시설)"),
    dict(name="vil", prefix="vil__", files=["마을체육시설_시군구_연도.csv"], sido="sido_clean", sgg="sigungu_clean", merge_counts=True,
         ratio=lambda c: False, label="문체부 전국공공체육시설현황(마을체육시설)"),
    dict(name="canc", prefix="canc__", files=["암발생_시군구_시기_wide.csv"], sido="sido", sgg="sigungu", merge_counts=True,
         ratio=lambda c: c.endswith("발생률"), label="암등록통계 시군구 암발생(5년 시기, year=종료연도)", optional=True),
    dict(name="eloc", prefix="eloc__", files=["F_e지방지표_시군구_연도_wide.csv"], sido="sido", sgg="sigungu", merge_counts=False,
         ratio=lambda c: True, label="e-지방지표", optional=True, pool="first", subs_unmatched=True, merged_unmatched=True),
    dict(name="etc", prefix="etc__", files=["E_기타_시군구_연도_wide.csv"], sido="sido", sgg="sigungu", merge_counts=False,
         ratio=lambda c: True, label="기타 시군구 공변량", optional=True, pool="first", subs_unmatched=True, merged_unmatched=True),
]


def load_chk_formulas():
    f = os.path.join(PAN, "NHIS검진_변수사전.csv")
    out = {}
    if not os.path.exists(f): return out
    d = pd.read_csv(f, dtype=str)
    for r in d[d.kind == "파생"].itertuples():
        m = re.match(r"^\s*([^/]+?)\s*/\s*([^*]+?)\s*\*\s*100\s*$", str(r.measure))
        if m: out[r.column] = (m.group(1).strip(), m.group(2).strip())
    for r in d[(d.unit == "%") & (d.kind != "파생")].itertuples():
        m = re.match(r"^(대사_.+)_비율$", str(r.column))
        if m: out[r.column] = (m.group(1) + "_인원", "대사_수검자수")
    return out


def load_source(src, master_keys, log):
    frames = []
    for f in src["files"]:
        p = os.path.join(PAN, f)
        if not os.path.exists(p):
            log["missing_files"].append(f); continue
        d = pd.read_csv(p, dtype={src["sido"]: str, src["sgg"]: str, "sido": str, "sigungu": str, "level": str, "sgg_code": str,
                                  "sgg_code_std": str, "parent_code": str, "code_note": str, "region_nm": str, "kchs_code": str,
                                  "sido_clean": str, "sigungu_clean": str, "sido_short": str, "sheet": str, "source_file": str}, low_memory=False)
        d["year"] = pd.to_numeric(d["year"], errors="coerce").astype("Int64")
        if "level" in d.columns and src.get("level"):
            n0 = len(d); d = d[d["level"] == src["level"]]; log["level_filter"][f] = (n0, len(d))
        d = d.rename(columns={src["sido"]: "_sido", src["sgg"]: "_sgg"}) if src["sido"] != "sido" or src["sgg"] != "sigungu" else d.rename(columns={"sido": "_sido", "sigungu": "_sgg"})
        d = d.drop(columns=[c for c in ID_DROP + ["sido", "sigungu"] if c in d.columns and c not in ("_sido", "_sgg")])
        frames.append(d)
    if not frames: return None, None
    d = frames[0]
    for f2 in frames[1:]:
        d = d.merge(f2, on=["year", "_sido", "_sgg"], how="outer", suffixes=("", "_dup"))
        d = d.drop(columns=[c for c in d.columns if c.endswith("_dup")])
    d = d[d["_sgg"].notna() & d["_sido"].notna()].copy()
    var_cols = [c for c in d.columns if c not in ("year", "_sido", "_sgg")]
    for c in var_cols: d[c] = pd.to_numeric(d[c], errors="coerce")

    pairs = d[["_sido", "_sgg"]].drop_duplicates()
    mp = {(a, b): map_region(a, b, master_keys, subs_unmatched=src.get("subs_unmatched", False), merged_unmatched=src.get("merged_unmatched", False)) for a, b in zip(pairs._sido, pairs._sgg)}
    d["key"] = [mp[(a, b)]["key"] for a, b in zip(d._sido, d._sgg)]
    d["_role"] = [mp[(a, b)]["role"] for a, b in zip(d._sido, d._sgg)]
    d["_flag"] = [mp[(a, b)]["flag"] for a, b in zip(d._sido, d._sgg)]
    d["_rmode"] = [mp[(a, b)]["ratio_mode"] or "" for a, b in zip(d._sido, d._sgg)]
    d["_nonnull"] = d[var_cols].notna().any(axis=1)

    um = d[d.key.isna()]
    um_tbl = um.groupby(["_sido", "_sgg"]).agg(years=("year", lambda s: f"{s.min()}~{s.max()}"), n_rows=("year", "size"),
                                                n_rows_with_values=("_nonnull", "sum")).reset_index()
    um_tbl.insert(0, "source", src["name"])
    log["unmatched"].append(um_tbl)
    d = d[d.key.notna()].copy()

    sp = d[d._role.ne("primary") | d._flag.ne("")].groupby(["_sido", "_sgg", "key", "_role", "_flag"]).agg(
        years=("year", lambda s: f"{s.min()}~{s.max()}"), n_rows=("year", "size")).reset_index()
    sp = sp.rename(columns={"_sido": "sido", "_sgg": "sigungu", "_role": "role", "_flag": "flag"})
    sp.insert(0, "source", src["name"]); log["special"].append(sp)

    ratio_fn = src["ratio"]
    formulas = load_chk_formulas() if ratio_fn == "dict" else {}
    if ratio_fn == "dict":
        chk_dict = pd.read_csv(os.path.join(PAN, "NHIS검진_변수사전.csv"), dtype=str).drop_duplicates("column").set_index("column")
        is_ratio = lambda c: (c in formulas) or (chk_dict.unit.get(c, "") == "%")
    else:
        is_ratio = ratio_fn
    count_cols = [c for c in var_cols if not is_ratio(c)] if src["merge_counts"] else []
    ratio_cols = [c for c in var_cols if c not in count_cols]

    g = d.groupby(["key", "year"], sort=False)
    sizes = g.size()
    multi_idx = sizes[sizes > 1].index
    single = d.set_index(["key", "year"]).loc[~d.set_index(["key", "year"]).index.isin(multi_idx)]
    single_out = single[var_cols].copy(); single_out["boundary_flag"] = single["_flag"].values
    single_out = single_out.reset_index()

    rows = []; overlap_notes = []
    if len(multi_idx):
        multi = d.set_index(["key", "year"]).loc[multi_idx]
        for (key, yr), grp in multi.groupby(level=[0, 1], sort=False):
            flags = sorted(set(f for f in grp._flag if f))
            roles = grp._role.tolist(); rmodes = set(m for m in grp._rmode if m)
            out = {}
            prim = grp[grp._role == "primary"]
            if len(prim) and (grp._role == "sub").any():
                grp = grp[grp._role != "sub"]; rmodes = set(m for m in grp._rmode if m)
                flags = sorted((set(flags) - {"하위구합산"}) | {"하위구제외"})
                if len(grp) == 1:
                    out = {c: grp[c].iloc[0] for c in var_cols}
                    out["key"] = key; out["year"] = yr; out["boundary_flag"] = ";".join(flags); rows.append(out); continue
            for c in count_cols:
                out[c] = grp[c].sum(min_count=1)
            if src.get("pool") == "first":
                ordered = pd.concat([prim, grp[grp._role != "primary"]])
                for c in var_cols:
                    v = ordered[c].dropna(); out[c] = v.iloc[0] if len(v) else np.nan
                out["key"] = key; out["year"] = yr; out["boundary_flag"] = ";".join(flags); rows.append(out); continue
            if src.get("pool") == "kchs" and len(prim) >= 2:
                pooled, did_pool = pool_kchs(prim, var_cols)
                out.update(pooled)
                if did_pool: flags = sorted(set(flags) | {"보건소합산"})
                out["key"] = key; out["year"] = yr; out["boundary_flag"] = ";".join(flags); rows.append(out); continue
            for c in ratio_cols:
                v = np.nan
                if c in formulas and all(x in out for x in formulas[c]):
                    num, den = formulas[c]
                    v = out[num] / out[den] * 100 if (pd.notna(out[den]) and out[den] != 0 and pd.notna(out[num])) else np.nan
                elif "nan" in rmodes:
                    v = np.nan
                else:
                    pv = prim[c].dropna() if len(prim) else pd.Series(dtype=float)
                    if len(pv) == 1: v = pv.iloc[0]
                    elif len(pv) == 0 and src["merge_counts"] is False:
                        nv = grp[c].dropna(); v = nv.iloc[0] if len(nv) == 1 else np.nan
                    else: v = np.nan
                out[c] = v
            if (np.array(roles) == "primary").sum() >= 2:
                ov = int((prim[var_cols].notna().sum(axis=0) >= 2).sum())
                overlap_notes.append(dict(source=src["name"], key=key, year=int(yr), rows="+".join(f"{a}/{b}" for a, b in zip(grp._sido, grp._sgg)),
                                          n_cols_overlap=ov))
            out["key"] = key; out["year"] = yr; out["boundary_flag"] = ";".join(flags)
            rows.append(out)
    agg = pd.concat([single_out, pd.DataFrame(rows)], ignore_index=True) if rows else single_out
    agg["year"] = agg["year"].astype(int)
    log["overlap"].extend(overlap_notes)
    meta = dict(n_rows_in=int(len(d)), var_cols=var_cols, count_cols=count_cols, ratio_cols=ratio_cols, n_multi=int(len(multi_idx)),
                flag_counts=agg["boundary_flag"].replace("", np.nan).dropna().str.split(";").explode().value_counts().to_dict(),
                match_by_year=agg.groupby("year").key.nunique().to_dict())
    return agg, meta


def pool_kchs(prim, var_cols):
    out = {}
    overlap = (prim[var_cols].notna().sum(axis=0) >= 2)
    if not overlap.any():
        for c in var_cols:
            v = prim[c].dropna(); out[c] = v.iloc[0] if len(v) else np.nan
        return out, False
    bases = sorted(set(c.rsplit("__", 1)[0] for c in var_cols if "__" in c))
    for b in bases:
        cN = f"{b}__N"; N = prim[cN] if cN in prim.columns else pd.Series(np.nan, index=prim.index)
        for suf in ["CR", "SR", "CR_SE", "SR_SE", "N"]:
            c = f"{b}__{suf}"
            if c not in prim.columns: continue
            v = prim[c]
            if v.notna().sum() <= 1:
                out[c] = v.dropna().iloc[0] if v.notna().any() else np.nan; continue
            if suf == "N": out[c] = v.sum(); continue
            m = v.notna() & N.notna()
            if m.sum() < 2 or N[m].sum() == 0: out[c] = np.nan; continue
            w = N[m] / N[m].sum()
            out[c] = float((w * v[m]).sum()) if not suf.endswith("SE") else float(np.sqrt((w ** 2 * v[m] ** 2).sum()))
    for c in var_cols:
        if c not in out:
            v = prim[c].dropna(); out[c] = v.iloc[0] if len(v) == 1 else np.nan
    return out, True


def dict_lookup(src_name):
    out = {}
    if src_name == "kchs" and _exists("KCHS_변수사전.csv"):
        d = pd.read_csv(os.path.join(PAN, "KCHS_변수사전.csv"), dtype=str)
        for r in d.itertuples(): out[r.column] = (r.TBL_ID, r.TBL_NM, r.ITM_NM, "명" if r.ITM_ID == "N" else "%")
    elif src_name == "nhis" and _exists("NHIS_변수사전.csv"):
        d = pd.read_csv(os.path.join(PAN, "NHIS_변수사전.csv"), dtype=str); d = d[d.output.str.contains("wide", na=False)]
        for r in d.itertuples(): out[r.column] = (r.TBL_ID, r.TBL_NM, " ".join(x for x in [r.C2_NM, r.ITM_NM] if isinstance(x, str)), r.UNIT_NM)
    elif src_name == "chk" and _exists("NHIS검진_변수사전.csv"):
        d = pd.read_csv(os.path.join(PAN, "NHIS검진_변수사전.csv"), dtype=str).drop_duplicates("column")
        for r in d.itertuples(): out[r.column] = (r.tbl_id, f"{r.topic}/{r.category}", r.measure, r.unit)
    elif src_name in ("mort", "pop") and _exists("사망_변수사전.csv"):
        d = pd.read_csv(os.path.join(PAN, "사망_변수사전.csv"), dtype=str)
        fn = "사망원인_시군구_연도_wide.csv" if src_name == "mort" else "인구동향_시군구_연도_wide.csv"
        for r in d[d.파일 == fn].itertuples(): out[r.column] = (r.TBL_ID, r.TBL_NM, f"{r.ITM_NM} {r.설명}".strip(), r.단위)
    elif src_name == "canc" and _exists("암발생_변수사전.csv"):
        d = pd.read_csv(os.path.join(PAN, "암발생_변수사전.csv"), dtype=str)
        for r in d.itertuples(): out[r.column] = (r.TBL_ID, r.TBL_NM, getattr(r, "note", ""), "명" if r.column.endswith("발생자수") else "십만명당")
    elif src_name in ("eloc", "etc"):
        fn = {"eloc": "F_e지방지표_변수사전.csv", "etc": "E_기타_변수사전.csv"}[src_name]
        if _exists(fn):
            d = pd.read_csv(os.path.join(PAN, fn), dtype=str)
            g = lambda r, *names: next((getattr(r, n) for n in names if hasattr(r, n) and isinstance(getattr(r, n), str)), "")
            for r in d.itertuples():
                out[g(r, "column")] = (g(r, "TBL_ID"), g(r, "TBL_NM"), g(r, "ITM_NM"), g(r, "UNIT", "UNIT_NM", "unit"))
    return out


KCHS_DISEASE_TBL = {"DT_DIA_DOCTOR", "DT_HYPER_DOCTOR", "DT_117075_H_MENTAL_DEPRESS_SYM", "DT_H_MENTAL_DEPRESS", "DT_117075_MEA_OBE", "DT_H_OBE_OBE", "DT_H_OR_INCONV"}
CHK_DEP_PAT = re.compile(r"^(판정_유질환자|판정_질환의심|판정_일반질환의심|판정_고혈압당뇨질환의심|유질환자_|질환의심_|대사_(?!수검자수)|인지기능)")


def is_dependent(src_name, col, tbl_id):
    if src_name in ("mort", "canc"): return "Y"
    if src_name == "nhis": return "Y" if "주요질환" in str(tbl_id) or re.match(r"^(고혈압|당뇨|치주질환|관절염|감염병|간질환|정신질환|.*암)__", col) else "N"
    if src_name == "chk": return "Y" if CHK_DEP_PAT.match(col) else "N"
    if src_name == "kchs":
        return "Y" if (str(tbl_id) in KCHS_DISEASE_TBL and (col.endswith("__CR") or col.endswith("__SR"))) else "N"
    return "N"


def correct_changwon_population(master):
    # KOSIS DT_1YL20651E, T20, 48120, 2023. The historical FBD input repeats the 2022 value.
    rows = master.index[master.key.eq("경상남도|창원시") & master.year.eq(2023)]
    if not len(rows) and not master.year.eq(2023).any():
        return
    if len(rows) != 1:
        raise ValueError("Missing or duplicated Changwon 2023")
    row = rows[0]
    if master.at[row, "population"] not in (1021487, 1009038):
        raise ValueError("Reconcile the new Changwon population vintage before use")
    field = "eloc__DT_1YL20651E__계"
    if field in master and pd.notna(master.at[row, field]) and float(master.at[row, field]) != 1009038:
        raise ValueError("Changwon population disagrees with the independent KOSIS total")
    master.at[row, "population"] = 1009038


def main():
    mk = pd.read_csv(os.path.join(PAN, "시군구_마스터키.csv"), dtype=str)
    mk = mk[mk.fbd_sido.notna()].drop_duplicates("key")
    master_keys = set(mk.key)
    assert len(master_keys) == 229, f"마스터키 229개 아님: {len(master_keys)}"
    keyinfo = mk.set_index("key")[["fbd_region"]]

    log = dict(missing_files=[], level_filter={}, unmatched=[], special=[], overlap=[], sources={})
    aggs = {}
    for src in SOURCES:
        if all(not _exists(f) for f in src["files"]):
            if not src.get("optional"): log["missing_files"] += src["files"]
            print(f"[skip] {src['name']}: 파일 없음"); continue
        agg, meta = load_source(src, master_keys, log)
        if agg is None: continue
        agg = agg.rename(columns={c: src["prefix"] + c for c in meta["var_cols"]})
        agg = agg.rename(columns={"boundary_flag": f"_flag_{src['name']}"})
        aggs[src["name"]] = (agg, meta, src)
        print(f"[ok] {src['name']}: 입력행 {meta['n_rows_in']} → 집계행 {len(agg)}, 변수 {len(meta['var_cols'])}, 합산그룹 {meta['n_multi']}")

    pop = pd.read_csv(FBD, usecols=["year", "sido", "sigungu", "population"], dtype={"sido": str, "sigungu": str}).drop_duplicates()
    pop["key"] = [make_key(a, b) for a, b in zip(pop.sido, pop.sigungu)]
    pop = pop.groupby(["key", "year"], as_index=False).population.first()
    pop_un = sorted(set(pop.key) - master_keys)

    years = sorted(set().union(*[set(a.year.unique()) for a, _, _ in aggs.values()]) | set(pop.year.unique()))
    grid = pd.MultiIndex.from_product([years, sorted(master_keys)], names=["year", "key"]).to_frame(index=False)
    grid["sido_std"] = grid.key.str.split("|").str[0]; grid["sigungu"] = grid.key.str.split("|").str[1]
    grid["fbd_region"] = grid.key.map(keyinfo.fbd_region)
    master = grid.merge(pop[["key", "year", "population"]], on=["key", "year"], how="left")
    flag_cols = []
    for name, (agg, meta, src) in aggs.items():
        master = master.merge(agg, on=["key", "year"], how="left"); flag_cols.append(f"_flag_{name}")
    fl = master[flag_cols].fillna("")
    master["boundary_flag"] = fl.apply(lambda r: ";".join(sorted(set(x for v in r for x in v.split(";") if x))), axis=1)
    master = master.drop(columns=flag_cols)
    front = ["year", "key", "sido_std", "sigungu", "fbd_region", "population", "boundary_flag"]
    master = master[front + [c for c in master.columns if c not in front]]
    master = master.sort_values(["year", "key"]).reset_index(drop=True)
    assert len(master) == 229 * len(years) and master.key.nunique() == 229 and not master.duplicated(["year", "key"]).any(), \
        f"마스터 격자 깨짐: {len(master)}행, 229×{len(years)} 기대"
    correct_changwon_population(master)
    master.to_parquet(os.path.join(PAN, "MASTER_시군구_연도_wide.parquet"), index=False)
    master.to_csv(os.path.join(PAN, "MASTER_시군구_연도_wide.csv.gz"), index=False, encoding="utf-8-sig", compression="gzip")

    rows = []
    for name, (agg, meta, src) in aggs.items():
        lk = dict_lookup(name)
        for c in meta["var_cols"]:
            mc = src["prefix"] + c; s = master[mc]
            tbl_id, tbl_nm, item, unit = lk.get(c, ("", src["label"], "", ""))
            if name in ("fac", "vil") and not unit: unit = "m2" if "면적" in c else "개소"
            nn = s.notna()
            yrs = master.loc[nn, "year"]
            y0, y1 = (int(yrs.min()), int(yrs.max())) if nn.any() else (None, None)
            n_last = int(master.loc[nn & (master.year == y1), "key"].nunique()) if nn.any() else 0
            rows.append(dict(column=mc, source=name, source_label=src["label"], orig_column=c, TBL_ID=tbl_id, TBL_NM=tbl_nm, item=item, unit=unit,
                             kind="합산가능" if c in meta["count_cols"] else ("비율·지표" if src["merge_counts"] else "비합산(표본조사·지표)"),
                             year_min=y0, year_max=y1, n_sgg_latest=n_last, n_nonnull=int(nn.sum()), dependent=is_dependent(name, c, tbl_id)))
    vd = pd.DataFrame(rows)
    vd.to_csv(os.path.join(PAN, "MASTER_변수사전.csv"), index=False, encoding="utf-8-sig")

    um = pd.concat(log["unmatched"], ignore_index=True) if log["unmatched"] else pd.DataFrame(columns=["source", "_sido", "_sgg", "years", "n_rows", "n_rows_with_values"])
    um = um.rename(columns={"_sido": "sido", "_sgg": "sigungu"})
    um.to_csv(os.path.join(PAN, "MASTER_병합_unmatched.csv"), index=False, encoding="utf-8-sig")
    sp = pd.concat(log["special"], ignore_index=True) if log["special"] else pd.DataFrame()

    nY = int((vd.dependent == "Y").sum()); nN = int((vd.dependent == "N").sum())
    L = [f"# MASTER 병합 점검 ({datetime.date.today()})", "",
         f"- 마스터: {len(master):,}행 × {master.shape[1]:,}열, 시군구 {master.key.nunique()}개 × 연도 {years[0]}~{years[-1]}({len(years)}개), 완전 격자",
         f"- 인구 분모(population): FBD 식중독최종.csv {int(pop.year.min())}~{int(pop.year.max())}, 마스터에 값 있는 행 {int(master.population.notna().sum()):,}" + (f" (FBD 키 미매칭: {pop_un})" if pop_un else ""),
         f"- 변수 {len(vd)}개 = 종속변수 후보(Y) {nY} / 공변량(N) {nN}", ""]
    if log["missing_files"]: L += ["## 없는 입력 파일(건너뜀)"] + [f"- {f}" for f in log["missing_files"]] + [""]
    L += ["## 원천별 열 수·집계 요약", "", "| 원천 | 접두어 | 입력 파일 | 변수 수 | 합산가능 열 | 비율·지표 열 | 집계행 | 합산 그룹(key×year) |", "|---|---|---|---|---|---|---|---|"]
    for name, (agg, meta, src) in aggs.items():
        L.append(f"| {src['label']} | {src['prefix']} | {', '.join(f for f in src['files'] if _exists(f))} | {len(meta['var_cols'])} | {len(meta['count_cols'])} | {len(meta['ratio_cols'])} | {len(agg):,} | {meta['n_multi']} |")
    L += ["", "## 원천별 매칭된 시군구 수(연도별)", ""]
    cov = pd.DataFrame({name: pd.Series(meta["match_by_year"]) for name, (agg, meta, src) in aggs.items()}).reindex(years).fillna(0).astype(int)
    L += ["| year | " + " | ".join(cov.columns) + " |", "|---|" + "---|" * len(cov.columns)]
    L += [f"| {y} | " + " | ".join(str(v) for v in cov.loc[y]) + " |" for y in cov.index]
    L += ["", "## 연도 커버리지(원천×연도, 해당 원천 변수 중 하나라도 값 있는 시군구 수)", ""]
    cov2 = {}
    for name, (agg, meta, src) in aggs.items():
        cols = [src["prefix"] + c for c in meta["var_cols"]]
        cov2[name] = master[master[cols].notna().any(axis=1)].groupby("year").key.nunique()
    cov2["population"] = master[master.population.notna()].groupby("year").key.nunique()
    cov2 = pd.DataFrame(cov2).reindex(years).fillna(0).astype(int)
    L += ["| year | " + " | ".join(cov2.columns) + " |", "|---|" + "---|" * len(cov2.columns)]
    L += [f"| {y} | " + " | ".join(str(v) for v in cov2.loc[y]) + " |" for y in cov2.index]
    L += ["", "## 미매칭 (sido, sigungu) — 버리지 않고 MASTER_병합_unmatched.csv 에도 저장(마스터에는 들어가지 않음)", ""]
    if len(um):
        L += ["원천별 미매칭 행 수: " + ", ".join(f"{k} {int(v)}행({int(um[um.source==k].n_rows_with_values.sum())}행 값 있음)" for k, v in um.groupby("source").n_rows.sum().items()), ""]
    if len(um):
        L += ["| 원천 | sido | sigungu | 연도 | 행 수 | 값 있는 행 |", "|---|---|---|---|---|---|"]
        L += [f"| {r.source} | {r.sido} | {r.sigungu} | {r.years} | {r.n_rows} | {r.n_rows_with_values} |" for r in um.itertuples()]
    else: L.append("- 없음")
    L += ["", "## 합산·개명·표기정정 처리한 (sido, sigungu) 목록", "", "| 원천 | sido | sigungu | → key | 역할 | 플래그 | 연도 | 행 수 |", "|---|---|---|---|---|---|---|---|"]
    L += [f"| {r.source} | {r.sido} | {r.sigungu} | {r.key} | {r.role} | {r.flag} | {r.years} | {r.n_rows} |" for r in sp.itertuples()] if len(sp) else ["| - | | | | | | | |"]
    L += ["", "## 경계 규칙 적용 건수(마스터 행 기준, boundary_flag)", ""]
    fc = master.boundary_flag.replace("", np.nan).dropna().str.split(";").explode().value_counts()
    L += [f"- {k}: {v}행" for k, v in fc.items()] or ["- 없음"]
    L += ["", "원천별 플래그 건수(집계행 기준):"] + [f"- {name}: {meta['flag_counts']}" for name, (agg, meta, src) in aggs.items()]
    L += ["", "## 합산 때문에 결측 처리한 비율 셀 수(합산 그룹 행 × 비율 열, 원천별)"]
    for name, (agg, meta, src) in aggs.items():
        if not meta["ratio_cols"] or not src["merge_counts"]: continue
        fl = agg[f"_flag_{name}"].fillna("")
        rc_ = [src["prefix"] + c for c in meta["ratio_cols"]]
        for f_ in ["창원통합", "하위구합산"]:
            sub_ = agg[fl.str.contains(f_)]
            if len(sub_): L.append(f"- {name} / {f_}: {len(sub_)}행 × 비율 {len(rc_)}열, 결측 셀 {int(sub_[rc_].isna().sum().sum())}")
    if log["level_filter"]:
        L += ["", "## level 필터(시군구만 사용)"] + [f"- {f}: {a}행 → {b}행" for f, (a, b) in log["level_filter"].items()]
    L += ["", "## 같은 key×year에 원천 행이 둘 이상이면서 둘 다 '정식' 행인 경우(옛 이름·새 이름 공존 등) — 값이 겹친 열 수"]
    ov = pd.DataFrame(log["overlap"])
    if len(ov):
        ovs = ov.groupby(["source", "rows"]).agg(years=("year", lambda s: f"{s.min()}~{s.max()}"), n=("year", "size"), max_overlap_cols=("n_cols_overlap", "max"), sum_overlap=("n_cols_overlap", "sum")).reset_index()
        L += ["", "| 원천 | 행 조합 | 연도 | 건수 | 겹친 열 최대 | 겹친 열 합 |", "|---|---|---|---|---|---|"]
        L += [f"| {r.source} | {r.rows} | {r.years} | {r.n} | {r.max_overlap_cols} | {r.sum_overlap} |" for r in ovs.itertuples()]
        L += ["", "겹친 열 합이 0이면 값이 서로 보완적(한쪽만 값 있음)이라 합산해도 이중계상 없음. 0보다 크면 해당 열은 두 행 값을 더한 것이므로 확인 필요."]
    else: L.append("- 없음")
    L += ["", "## 처리 규칙 요약",
          "- 지역 키: 40_region_crosswalk.make_key(군위군→대구시, 연기군→세종시, 인천 남구→미추홀구, 제주/서귀포) 뒤 RENAME(1:1 개명·표기 정정) → MERGE(폐지 지역 합산) → 하위 구 규칙(고양·성남·수원·안산·안양·용인·부천·청주·천안·전주·포항·창원·화성의 구).",
          "- 합산 가능한 양(인원·일수·금액·건수·시설 수·면적·사망자 수·출생아 수)은 합, 비율·율·평균은 (1) 검진 파생비율은 합산 뒤 공식으로 재계산, (2) 창원통합은 결측, (3) 청원군·북제주군·남제주군 합산은 시 값 유지, (4) 하위구 합산은 결측.",
          "- KCHS는 표본조사 추정치라 합산하지 않음: 청원군(2008~2013) 행은 쓰지 않고 청주시 값 유지, boundary_flag='청원군합산' 기록. 같은 시군구에 보건소 코드 행이 둘(제주시·서귀포시)이면 응답자수(N) 가중 평균으로 풀링(N 합, SE는 가중 제곱합의 제곱근), boundary_flag='보건소합산'. 시도명 변형 행(강원도/강원특별자치도 등)은 열이 보완적이라 그대로 합침.",
          "- 같은 key×year에 시 행과 하위 구 행이 함께 있으면(암발생 등) 시 행만 쓰고 하위 구는 버림(boundary_flag='하위구제외'), 비율은 시 값 유지.",
          "- 암발생(canc__)은 5년 시기 종료연도(2003·2008·2013·2018) 행에만 값이 있고 다른 연도는 결측. 연기군(~2013 시기)은 crosswalk대로 세종시 키에 들어감. 청원군(~2013)·남북제주군(~2008)은 발생자수 합산·발생률은 시 값, 마산시·진해시(~2013)는 창원시에 발생자수 합산·발생률 결측.",
          "- e-지방지표(eloc__)는 같은 시군구가 시도명 변형(경기/경기도, 강원/강원도/강원특별자치도 등)으로 여러 행에 나뉘어 있어 make_key로 표준화한 뒤 같은 key×year를 열별 첫 비결측(정식 행 우선)으로 합침. 통합시 하위 구·소계·출장소·'(변동전)' 꼬리표 행은 합산하지 않고 미매칭으로만 보냄(단 '(변동전)' 꼬리표는 떼고 1:1 개명 지역 여주군·당진군·연기군→세종시·광주군·양주군·포천군·화성군·군위군·인천 남구·창원시·통합청주시는 기존 경계 규칙으로 매핑; 합산이 필요한 청원군·마산시·진해시·북제주군·남제주군은 지표를 합산할 수 없어 미매칭으로 보냄 — 청원군 값을 청주시로 쓰면 2004~2014에 군 단위 값 336셀이 섞여서 제외). '통합청주시'→청주시. 연도 범위(1992~2042, 추계인구 포함)가 마스터 격자의 연도 수를 늘림.",
          "- 사망원인·인구동향은 level=='시군구' 행만 사용(하위구 제외), 2010~ '통합창원시'→창원시, 2014~ 청주시는 원천이 이미 통합값.",
          "- 같은 key×year에 옛 이름·새 이름 행이 공존하면(여주군/여주시 등) 열별로 값이 한쪽에만 있는지 위 표로 점검.",
          "- 종속변수 후보(Y): 사망 전부, NHIS 주요질환, 검진 판정·유질환자·질환의심·대사증후군·인지기능, KCHS 당뇨·고혈압 진단경험률·우울·비만·저작불편(CR/SR). 그 외는 N."]
    open(os.path.join(PAN, "MASTER_병합_점검.md"), "w", encoding="utf-8").write("\n".join(L))

    print(f"\n마스터 {len(master):,}행 × {master.shape[1]:,}열 | 종속변수 후보(Y) {nY} / 공변량(N) {nN}")
    print("미매칭:", len(um), "건 →", um[["source", "sido", "sigungu"]].astype(str).agg("/".join, axis=1).tolist())
    return master, vd


if __name__ == "__main__":
    main()
