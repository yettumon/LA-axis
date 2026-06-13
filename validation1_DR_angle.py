"""
검증 1: Dipole Repeller — LA-axis와 각거리 계산
=================================================
논문: LA-axis 도출 — GCOD 주축 독립 논문
섹션: 4.1절 | Kim, Gun-Sik (Yettumon) 2026

핵심 결과:
  DR(RA=337.5°, Dec=0°) vs LA-axis(RA=330°, Dec=0°) = 7.5°
  우연 일치 확률 (random alignment probability) = 0.0043

주의사항:
  - 이 확률은 통계 검정 p-value가 아닌 random alignment probability임
  - DR 방향 중심값 사용 (Hoffman et al. 2017); 오차 전파는 논문 본문 참조
  - "물리적 연관성"이 아닌 "비무작위 정렬"을 의미함

데이터: 불필요 (순수 계산)
"""

import numpy as np

print("=" * 55)
print("검증 1: Dipole Repeller — LA-axis 각거리")
print("논문 4.1절 | Kim, Gun-Sik (Yettumon) 2026")
print("=" * 55)

# ── 좌표 정의 ──
LA_ra,  LA_dec = 330.0, 0.0   # LA-axis (RA=330°, Dec=0°)
DR_ra,  DR_dec = 337.5, 0.0   # Dipole Repeller 중심값 (Hoffman et al. 2017)
# 주의: DR 방향에는 관측 오차가 존재함.
# 논문 본문: "using the central value of the Dipole Repeller direction"
# 오차 전파 분석은 별도 sensitivity analysis 참조.

def angular_distance(ra1, dec1, ra2, dec2):
    """구면 코사인 법칙으로 각거리 계산 (degrees)"""
    r1, d1 = np.radians(ra1), np.radians(dec1)
    r2, d2 = np.radians(ra2), np.radians(dec2)
    cos_theta = (np.sin(d1)*np.sin(d2) +
                 np.cos(d1)*np.cos(d2)*np.cos(r1 - r2))
    return np.degrees(np.arccos(np.clip(cos_theta, -1, 1)))

angle = angular_distance(LA_ra, LA_dec, DR_ra, DR_dec)

# ── 우연 일치 확률: 구면 균일 분포에서 θ 이내 확률 ──
# P(θ) = (1 - cos(θ)) / 2  (등방 분포 가정, 단측)
# 통계 검정 p-value가 아닌 random alignment probability임
prob_random = (1 - np.cos(np.radians(angle))) / 2

print(f"\n[입력값]")
print(f"  LA-axis:         RA={LA_ra}°, Dec={LA_dec}°")
print(f"  Dipole Repeller: RA={DR_ra}°, Dec={DR_dec}° (중심값)")

print(f"\n[결과]")
print(f"  각거리:               {angle:.2f}°")
print(f"  Random alignment prob: {prob_random:.4f}")

# 논문 기재값 검증
expected_angle = 7.5
expected_prob  = 0.0043
ok_angle = abs(angle - expected_angle) < 0.1
ok_prob  = abs(prob_random - expected_prob) < 0.001

print(f"\n[논문 기재값 검증]")
print(f"  각거리: {angle:.2f}° (논문: {expected_angle}°) {'✅' if ok_angle else '❌'}")
print(f"  확률:   {prob_random:.4f} (논문: {expected_prob}) {'✅' if ok_prob else '❌'}")

print(f"\n[해석]")
print(f"  등방 분포에서 두 방향이 {angle:.1f}° 이내로 우연히 일치할 확률: {prob_random:.4f}")
print(f"  → 등방 무작위 분포에서 발생하기 어려운 정렬")
print(f"     (statistically unlikely to arise from isotropic random orientations)")
print(f"  → Cosmicflows-2 위너필터 재구성과 방법론적으로 독립된 수렴")
print(f"  ※ 이 수치는 random alignment probability이며,")
print(f"     물리적 인과관계(causation)를 직접 증명하지는 않음")
