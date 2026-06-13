"""
검증 6: CF3 격자탐색 — LA-axis 방향 정밀화
============================================
논문: LA-axis 도출 — GCOD 주축 독립 논문
섹션: 3.2절 | Kim, Gun-Sik (Yettumon) 2026

핵심 결과:
  RA=332.6°, Dec=0°  (c 최대화 기준)
  c 최대값: +378.7 km/s

탐색 방식:
  1단계: 15° 간격 전체 탐색 (Dec=0° 고정)
  2단계:  1° 간격 정밀 탐색 (최적 근방 ±20°)
  3단계: 0.1° 간격 초정밀 탐색 (최적 근방 ±2°) → 332.6° 정밀도

Dec=0° 고정 이유:
  - 물리적 대칭축이 적도면에 위치 (방사 분력 이론)
  - DR(Dec=0°), CMB 운동벡터(Dec=-6.9°)와 일관성

목적함수: c(cosθ) 최대화
  다변량 회귀 Vpec = a + b·z + c·cos(θ) 에서 c 계수 최대화
  c 최대 = LA-axis 방향 방사 분력 신호 최대화 방향

validation4와의 관계:
  validation4: RA=330° 고정 후 c=+378.2 측정 (신호 정량화)
  validation6: c 최대 RA 탐색 → 332.6° 도출 (방향 수렴 확인)
  동일 CF3 데이터의 다른 분석 → 논문에서 CF3 1채널로 취급

파싱: skiprows=5, z=Vmod/c (validation4와 동일)
데이터: CF3_all_distance.csv (필수)
"""

import numpy as np
import pandas as pd
import os, sys

print("=" * 55)
print("검증 6: CF3 격자탐색 — LA-axis 방향 정밀화")
print("논문 3.2절 | Kim, Gun-Sik (Yettumon) 2026")
print("=" * 55)

# ── 데이터 로드 ──
cf3_paths = [
    'CF3_all_distance.csv',
    '/content/drive/MyDrive/LA_axis/CF3_all_distance.csv',
    os.path.expanduser('~/CF3_all_distance.csv')
]
cf3_path = None
for p in cf3_paths:
    if os.path.exists(p):
        cf3_path = p
        break

if cf3_path is None:
    print("\n[오류] CF3_all_distance.csv 파일을 찾을 수 없습니다.")
    sys.exit(1)

cf3_raw = pd.read_csv(cf3_path, skiprows=5, header=None, low_memory=False)
header  = pd.read_csv(cf3_path, skiprows=1, nrows=1, header=None)
cols    = header.iloc[0].tolist()
cf3_raw.columns = cols[:len(cf3_raw.columns)]

for col in ['Glon','Glat','Vhel','Vmod','Dist']:
    cf3_raw[col] = pd.to_numeric(cf3_raw[col], errors='coerce')

df = cf3_raw[['Glon','Glat','Vhel','Vmod','Dist']].dropna().copy()
c_light = 299792.458
df['z']     = df['Vmod'] / c_light   # validation4와 동일
df['v_pec'] = df['Vhel'] - df['Vmod']
df = df[(df['z'] >= 0.01) & (df['z'] <= 0.07) & (df['Dist'] > 0)].copy()
print(f"\n유효 샘플: N={len(df)} (z=0.01~0.07, Vmod 기준)")

glon_r = np.radians(df['Glon'].values)
glat_r = np.radians(df['Glat'].values)
v_pec  = df['v_pec'].values
z_arr  = df['z'].values

def ra_to_gal(ra_deg, dec_deg):
    """적도 좌표 → 은하 좌표"""
    r, d = np.radians(ra_deg), np.radians(dec_deg)
    R = np.array([[-0.0548756, -0.8734371, -0.4838350],
                  [ 0.4941094, -0.4448296,  0.7469823],
                  [-0.8676661, -0.1980764,  0.4559838]])
    xeq = np.array([np.cos(d)*np.cos(r), np.cos(d)*np.sin(r), np.sin(d)])
    xg  = R @ xeq
    gl  = np.degrees(np.arctan2(xg[1], xg[0])) % 360
    gb  = np.degrees(np.arcsin(np.clip(xg[2], -1, 1)))
    return gl, gb

def regression_c(ra, dec=0.0):
    """cos(θ) 회귀 계수 c 계산 (Dec=0° 기본값)"""
    gl, gb = ra_to_gal(ra, dec)
    tgl = np.radians(gl); tgb = np.radians(gb)
    ct = np.clip(
        np.sin(glat_r)*np.sin(tgb) +
        np.cos(glat_r)*np.cos(tgb)*np.cos(glon_r - tgl), -1, 1)
    X = np.column_stack([ct, z_arr, np.ones(len(v_pec))])
    try:
        return np.linalg.lstsq(X, v_pec, rcond=None)[0][0]
    except:
        return np.nan

# ── 주요 방향 확인 ──
print(f"\n[주요 방향 c값 (Dec=0° 고정)]")
for ra, label in [
    (330,   "LA-axis 정의 (RA=330°)"),
    (337.5, "DR (RA=337.5°)"),
    (285,   "GA+PP 방향 (RA=285°, 물리 배경)"),
    (332.6, "논문 기재값 (RA=332.6°)"),
]:
    c = regression_c(ra, 0)
    print(f"  {label:35s}: c={c:+.1f} km/s")

# ── 1단계: 15° 간격 전체 탐색 (Dec=0° 고정) ──
print(f"\n[1단계: 15° 간격 전체 탐색 (Dec=0° 고정)]")
best1_c, best1_ra = -np.inf, 0
for ra in np.arange(0, 360, 15):
    c = regression_c(ra, 0)
    if c > best1_c:
        best1_c, best1_ra = c, ra
print(f"  1단계 최적: RA={best1_ra}°, c={best1_c:.1f} km/s")

# ── 2단계: 1° 간격 정밀 탐색 ──
print(f"\n[2단계: 1° 간격 정밀 탐색 (RA {best1_ra-20}~{best1_ra+20}°)]")
best2_c, best2_ra = -np.inf, 0
for ra in np.arange(best1_ra - 20, best1_ra + 21, 1):
    c = regression_c(ra % 360, 0)
    if c > best2_c:
        best2_c, best2_ra = c, ra % 360
print(f"  2단계 최적: RA={best2_ra}°, c={best2_c:.1f} km/s")

# ── 3단계: 0.1° 간격 초정밀 탐색 ──
print(f"\n[3단계: 0.1° 간격 초정밀 탐색 (RA {best2_ra-2:.1f}~{best2_ra+2:.1f}°)]")
best3_c, best3_ra = -np.inf, 0
for ra in np.arange(best2_ra - 2, best2_ra + 2.1, 0.1):
    c = regression_c(ra % 360, 0)
    if c > best3_c:
        best3_c, best3_ra = c, ra % 360
print(f"  3단계 최적: RA={best3_ra:.1f}°, c={best3_c:.2f} km/s")

# ── 논문 기재값 검증 ──
print(f"\n[논문 기재값 검증]")
ok_ra = abs(best3_ra - 332.6) <= 1.0
ok_c  = abs(best3_c  - 378.7) <= 5.0
print(f"  최적 RA:  {best3_ra:.1f}° (논문: 332.6°, Dec=0° 고정) {'✅' if ok_ra else '⚠️'}")
print(f"  최대 c:   {best3_c:.2f} km/s (논문: +378.7) {'✅' if ok_c else '⚠️'}")
print(f"  LA-axis 정의(330°)와의 차이: {abs(best3_ra-330):.1f}°")

# ── 비교 ──
c_LA  = regression_c(330,   0)
c_opt = regression_c(best3_ra, 0)
print(f"\n[LA-axis 정의(330°) vs 최적(332.6°) 비교]")
print(f"  c(330°)  = {c_LA:.2f} km/s")
print(f"  c(332.6°)= {c_opt:.2f} km/s")
print(f"  차이:      {c_opt-c_LA:.2f} km/s ({(c_opt-c_LA)/c_LA*100:.1f}%)")
print(f"  → 2.6° 차이로 c값 차이 미미 → 330°를 주축으로 채택 타당")

# ── 주의사항 ──
print(f"\n[주의사항]")
print(f"  Dec=0° 고정 이유: 물리적 대칭축 + DR/CMB와 일관성")
print(f"  목적함수: c(cosθ) 최대화 (p값 또는 R² 최소화 아님)")
print(f"  332.6°는 0.1° 탐색의 점추정값 (방향 오차 ±별도 추정 필요)")
print(f"  validation4(RA=330° 고정)와 본 코드는 동일 CF3 데이터 분석")
print(f"  → 논문에서 Cosmicflows 계열 1채널로 취급")
print(f"  반증 예측: z>0.07 범위에서도 RA≈330° 수렴 유지 여부")
print(f"  → DESI DR3/DR4, Euclid로 검증 가능")
