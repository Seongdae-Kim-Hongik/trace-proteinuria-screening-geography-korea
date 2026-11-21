#!/bin/bash
# run_A_gate.sh - Assembles the analytic panel (80) and runs the optional data and panel checks (ST_GATE_DIR). No model fitting.
S="${ST_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"; cd "$S"; T=NBKID; D=261121
ST_VER=std3 python3 01_스크립트/80_build_disease_panel_v.py $T > 04_분석/$T/run_80.log 2>&1; echo "80 rc=$?"
[ -z "${ST_GATE_DIR:-}" ] || python3 "$ST_GATE_DIR/data_gate.py" 04_분석/$T/panel_${T}_std3.csv --region fbd_region > 04_분석/$T/data_gate_${T}_std3_$D.txt 2>&1; echo "data_gate rc=$?"
G="${ST_GATE_DIR:-}"
[ -z "$G" ] || python3 "$G/panel_gate.py" 04_분석/$T/panel_${T}_std3.csv --name-col fbd_region --cols population,E,cases,offset --positive population,E,offset --years 2011-2024 > 04_분석/$T/panel_gate_${T}_std3_$D.txt 2>&1; echo "panel_gate(panel) rc=$?"
[ -z "$G" ] || python3 "$G/panel_gate.py" 02_패널/NBKID_Estd3.csv --name-col fbd_region --cols examinees,cases,E_std,E_std_age,E_pop,E_crude --positive examinees,E_std,E_std_age,E_pop,E_crude --years 2011-2024 > 04_분석/$T/panel_gate_${T}_Estd3_$D.txt 2>&1; echo "panel_gate(E) rc=$?"
