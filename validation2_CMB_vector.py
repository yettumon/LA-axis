"""
검증 2: CMB 운동벡터 분석 (Planck 2018)
=========================================
논문: LA-axis 도출 — GCOD 주축 독립 논문
섹션: 4.2절 | Kim, Gun-Sik (Yettumon) 2026

분석 내용:
  1. CMB dipole 방향 → J2000 변환
  2. CMB 운동 방향(RA=167.9°)과 LA-axis Leo 끝(RA=150°) 각거리 (19.1°)
  3. CMB 운동벡터의 LA-axis 방향 투영 (94.5%)

핵심 결과:
  CMB 운동 방향(RA=167.9°, Dec=-6.9°)과 LA-axis Leo 끝(RA=150°, Dec=0°)
  각거리 = 19.1°, random alignment probability = 0.0275
  부호 있는 투영: dot = -0.945 (반평행 → CMB는 Leo 반대, 즉 Aquarius 방향)
  투영 크기: |dot| = 0.945 = 94.5%, 투영 속도 = 348.6 km/s

주의:
  CMB 결과는 Fisher 결합에 미포함.
  정성적 독립 확인(4.2절)으로만 활용.
  LA-axis는 물리 원리에서 선험적 정의(2.1절) — CMB 데이터와 무관.

데이터: 불필요 (Planck 2018 공개값 사용)
"""

import numpy as np

print("=" * 55)
print("검증 2: CMB 운동벡터 분석")
print("논문 4.2절 | Kim, Gun-Sik (Yettumon) 2026")
print("=" * 55)

# ── Planck 2018 CMB dipole (galactic 좌표) ──
# l=264.021°, b=+48.253° (Planck 2018 Table 1)
cmb_l = 264.021
cmb_b =  48.253
v_cmb = 369.0  # km/s (태양계 CMB 운동 속도)

# ── 은하 좌표 → J2000 적도 좌표 변환 ──
def gal_to_j2000(l_deg, b_deg):
    l, b = np.radians(l_deg), np.radians(b_deg)
    R = np.array([
        [-0.0548756, -0.8734371, -0.4838350],
        [ 0.4941094, -0.4448296,  0.7469823],
        [-0.8676661, -0.1980764,  0.4559838]
    ])
    x_gal = np.array([np.cos(b)*np.cos(l),
                       np.cos(b)*np.sin(l),
                       np.sin(b)])
    x_eq = R.T @ x_gal
    dec = np.degrees(np.arcsin(np.clip(x_eq[2], -1, 1)))
    ra  = np.degrees(np.arctan2(x_eq[1], x_eq[0])) % 360
    return ra, dec

cmb_ra, cmb_dec = gal_to_j2000(cmb_l, cmb_b)

print(f"\n[CMB dipole 방향 (Planck 2018)]")
print(f"  은하 좌표: l={cmb_l}°, b={cmb_b}°")
print(f"  J2000 변환 (운동 방향): RA={cmb_ra:.1f}°, Dec={cmb_dec:.1f}°")

# ── LA-axis Leo 끝 (RA=150°, Dec=0°) ──
# 본 분석은 CMB 운동 방향(RA=167.9°)과 Leo 끝(RA=150°) 직접 비교
LA_leo_ra  = 150.0
LA_leo_dec = 0.0

def angular_distance(ra1, dec1, ra2, dec2):
    r1, d1 = np.radians(ra1), np.radians(dec1)
    r2, d2 = np.radians(ra2), np.radians(dec2)
    cos_t = (np.sin(d1)*np.sin(d2) +
             np.cos(d1)*np.cos(d2)*np.cos(r1 - r2))
    return np.degrees(np.arccos(np.clip(cos_t, -1, 1)))

angle = angular_distance(cmb_ra, cmb_dec, LA_leo_ra, LA_leo_dec)

# ── 우연 일치 확률 (random alignment probability) ──
# 통계 검정 p-value가 아닌 구면 균일 분포 기반 기하학적 확률
prob_random = (1 - np.cos(np.radians(angle))) / 2

print(f"\n[CMB 운동 방향과 Leo 끝(RA=150°) 각거리]")
print(f"  CMB 운동 방향: RA={cmb_ra:.1f}°, Dec={cmb_dec:.1f}°")
print(f"  LA-axis Leo:   RA={LA_leo_ra}°, Dec={LA_leo_dec}°")
print(f"  각거리:               {angle:.1f}°")
print(f"  Random alignment prob: {prob_random:.4f}")

# ── CMB 벡터의 LA-axis 방향 투영 ──
# LA-axis 주방향 단위벡터 (RA=330°, Dec=0°, Aquarius 방향)
LA_main_ra = 330.0
la_vec = np.array([
    np.cos(np.radians(0.0))*np.cos(np.radians(LA_main_ra)),
    np.cos(np.radians(0.0))*np.sin(np.radians(LA_main_ra)),
    np.sin(np.radians(0.0))
])
# CMB 운동 방향 단위벡터
cmb_vec = np.array([
    np.cos(np.radians(cmb_dec))*np.cos(np.radians(cmb_ra)),
    np.cos(np.radians(cmb_dec))*np.sin(np.radians(cmb_ra)),
    np.sin(np.radians(cmb_dec))
])

# 부호 있는 투영 및 크기 모두 출력
dot_signed  = np.dot(la_vec, cmb_vec)
projection  = abs(dot_signed)
v_projected = v_cmb * projection

print(f"\n[CMB 운동벡터 LA-axis 방향 투영]")
print(f"  CMB 속도:       {v_cmb} km/s")
print(f"  부호 있는 투영: {dot_signed:.4f}")
print(f"  (음수 → CMB가 Aquarius 반대, 즉 Leo 방향 운동)")
print(f"  투영 크기:      |dot| = {projection:.4f} = {projection*100:.1f}%")
print(f"  투영 속도:      {v_projected:.1f} km/s")

# ── 논문 기재값 검증 ──
print(f"\n[논문 기재값 검증]")
print(f"  각거리: {angle:.1f}° (논문: 19.1°) {'✅' if abs(angle-19.1)<0.5 else '❌'}")
print(f"  확률:   {prob_random:.4f} (논문: 0.0275) {'✅' if abs(prob_random-0.0275)<0.002 else '❌'}  ※ random alignment probability")
print(f"  투영:   {projection*100:.1f}% (논문: 94.5%) {'✅' if abs(projection*100-94.5)<1 else '❌'}")
print(f"  투영속도: {v_projected:.1f} km/s (논문: 348.6 km/s) {'✅' if abs(v_projected-348.6)<5 else '❌'}")

print(f"\n[주의]")
print(f"  CMB 결과(p=0.0275)는 Fisher 결합에 미포함.")
print(f"  논문 4.2절의 정성적 독립 확인으로만 활용.")
print(f"  LA-axis는 CMB 데이터와 무관하게 물리 원리에서 선험적 정의(2.1절).")
