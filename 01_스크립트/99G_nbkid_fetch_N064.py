# 99G_nbkid_fetch_N064.py - Downloads KOSIS tables for the judgement check: DT_35007_N064 (urine dipstick protein result by province and sex, 2011-2024)
#  and DT_35007_N065 (by sex and age, 2011). The API key is read from the KOSIS_API_KEY environment variable.
#  Output: 00_원천/C_NHIS_검진보험/KOSIS_NBKID_261121/*.json
import os, json, urllib.request, urllib.parse
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(R, "00_원천", "C_NHIS_검진보험", "KOSIS_NBKID_261121"); os.makedirs(OUT, exist_ok=True)
K = os.environ["KOSIS_API_KEY"]
for t, yrs in [("DT_35007_N064", range(2011, 2025)), ("DT_35007_N065", range(2011, 2012))]:
    p = os.path.join(OUT, f"{t}_{min(yrs)}_{max(yrs)}.json")
    if os.path.exists(p): continue
    out = []
    for y in yrs:
        u = "https://kosis.kr/openapi/Param/statisticsParameterData.do?" + urllib.parse.urlencode(dict(method="getList", apiKey=K, orgId="350", tblId=t, itmId="ALL", objL1="ALL", objL2="ALL", objL3="ALL", format="json", jsonVD="Y", prdSe="Y", startPrdDe=str(y), endPrdDe=str(y)))
        d = json.loads(urllib.request.urlopen(u, timeout=120).read().decode()); out += d if isinstance(d, list) else []; print(t, y, len(d) if isinstance(d, list) else d)
    json.dump(out, open(p, "w"), ensure_ascii=False)
