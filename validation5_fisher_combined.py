"""
검증 5: Fisher 결합 유의도 - v23.2 최종
=========================================
논문: LA-axis 도출 - GCOD 주축 독립 논문
섹션: 6절 | Kim, Gun-Sik (Yettumon) 2026

채널 구성 (v23.2 최종 — 실질 독립 3채널):
  DR / CF3 / Pantheon+ = 3채널
  CMB는 Fisher 제외 (정성 확인만)
  채널 1: Dipole Repeller (Cosmicflows-2, 위너필터 재구성)
          p = 0.0043
  채널 2: CF3 cos(theta) 회귀 (Cosmicflows-3, TF/SBF)
          p = 3.31e-4
  채널 3: Pantheon+SH0ES 캘리브레이터 편중 (이항검정)
          보수적(귀무=50%):   p = 4.11e-5
          물리적(귀무=33.7%): p = 3.09e-12

  * CMB 운동벡터 (p=0.0275): Fisher 미포함, 정성 확인만

핵심 결과:
  보수적 (귀무=50%):   p = 1.77e-8  -> 5.5sigma
  물리적 (귀무=33.7%): p = 3.69e-15 -> 7.8sigma

독립성 근거:
  DR: Cosmicflows-2, 위너필터, 속도장
  CF3: Cosmicflows-3, TF/SBF, 특이속도 회귀
  Cal: Pantheon+, 이항검정, 거리지수
  -> 데이터셋, 방법론, 물리량 모두 상이

데이터: 불필요 (p값 직접 입력)
"""

import numpy as np
from scipy.stats import norm, chi2

print("=" * 55)
print("검증 5: Fisher 결합 유의도 (v23.2)")
print("논문 6절 | Kim, Gun-Sik (Yettumon) 2026")
print("=" * 55)

# ── 채널별 p값 ──
p_dr       = 0.0043   # 채널 1: DR (random alignment probability - 기하학적 확률)
p_cf3      = 3.31e-4  # 채널 2: CF3 cos(theta) 회귀 (통계 검정 p-value)
p_cal_cons = 4.11e-5  # 채널 3 보수적 (귀무=50%) (통계 검정 p-value)
p_cal_phys = 3.09e-12 # 채널 3 물리적 (귀무=33.7%) (통계 검정 p-value)
p_cmb      = 0.0275   # CMB (Fisher 미포함)

print(f"\n[채널별 p값]")
print(f"  채널 1 (DR, 7.5°):            p = {p_dr:.4f}  ({norm.isf(p_dr):.2f}sigma)")
print(f"  채널 2 (CF3 회귀, c=+378.2):  p = {p_cf3:.2e}  ({norm.isf(p_cf3):.2f}sigma)")
print(f"  채널 3 보수적 (귀무=50%):     p = {p_cal_cons:.2e}  ({norm.isf(p_cal_cons):.2f}sigma)")
print(f"  채널 3 물리적 (귀무=33.7%):   p = {p_cal_phys:.2e}  ({norm.isf(p_cal_phys):.2f}sigma)")
print(f"  CMB (미포함, 참고):            p = {p_cmb:.4f}  ({norm.isf(p_cmb):.2f}sigma)")

# 주의: 채널 1(DR)은 random alignment probability, 채널 2-3은 통계 검정 p-value
# Fisher 방법은 uniform[0,1] 분포 가정
# random alignment prob도 귀무가설 하에서 uniform[0,1] 충족 -> 혼합 가능
def fisher_combined(p_list, label):
    """Fisher 결합 통계: chi2 = -2 * sum(ln(p_i)), df=2k"""
    chi2_stat = -2 * np.sum(np.log(p_list))
    df = 2 * len(p_list)
    p_combined = chi2.sf(chi2_stat, df)
    sigma = norm.isf(p_combined)
    print(f"\n  [{label}]")
    print(f"    chi2 = {chi2_stat:.2f} (df={df})")
    print(f"    p    = {p_combined:.2e}")
    print(f"    sigma = {sigma:.1f}")
    return p_combined, sigma

print(f"\n[Fisher 결합 결과]")
p_cons, sig_cons = fisher_combined(
    [p_dr, p_cf3, p_cal_cons], "3채널 보수적 (귀무=50%)")
p_phys, sig_phys = fisher_combined(
    [p_dr, p_cf3, p_cal_phys], "3채널 물리적 (귀무=33.7%)")

# CMB 포함 참고용
print(f"\n  [참고: CMB 포함 4채널]")
fisher_combined([p_dr, p_cf3, p_cal_cons, p_cmb], "4채널 보수적")
fisher_combined([p_dr, p_cf3, p_cal_phys, p_cmb], "4채널 물리적")

# ── 논문 기재값 검증 ──
print(f"\n[논문 기재값 검증]")
ok_cons = abs(sig_cons - 5.5) < 0.1
ok_phys = abs(sig_phys - 7.8) < 0.1
print(f"  보수적: p={p_cons:.2e}, {sig_cons:.1f}sigma (논문: 1.77e-8, 5.5sigma) {'✅' if ok_cons else '❌'}")
print(f"  물리적: p={p_phys:.2e}, {sig_phys:.1f}sigma (논문: 3.69e-15, 7.8sigma) {'✅' if ok_phys else '❌'}")

# ── 구버전 비교 참고 ──
print(f"\n[참고: 채널 구성 변경 이력 비교]")
fisher_combined([p_dr, p_cal_cons], "v22 이전: DR+CF3 통합 2채널 보수적")
fisher_combined([p_dr, p_cal_phys], "v22 이전: DR+CF3 통합 2채널 물리적")

# ── 독립성 근거 ──
print(f"\n[독립성 근거 요약 (논문 표3, 6절 [124]절)]")
print(f"  채널 1 (DR):  Cosmicflows-2 / 위너필터 유체역학 / 속도장")
print(f"  채널 2 (CF3): Cosmicflows-3 / TF+SBF 직접 측정 / 특이속도 회귀")
print(f"  채널 3 (Cal): Pantheon+SH0ES / 이항검정 / 거리지수")
print(f"  -> 데이터셋, 방법론, 물리량 3중 독립")
print(f"  ※ RA=330 채택 근거: CMB dipole + Axis of Evil + Great Attractor 수렴 (논문 3.2절)")
print(f"     DR/CF3는 채택 이후 독립 검증 -- RA=330 선정에 미사용")
print(f"  ※ Fisher 방법: 귀무가설 하 p~uniform[0,1] 가정")
print(f"     random alignment prob(채널1)도 귀무가설 하 uniform[0,1] 충족")
print(f"     -> 이질적 p값 혼합에도 Fisher 적용 가능")

print(f"\n[채널 구성 변경 이력 - 심사자 방어]")
print(f"  구버전: DR+CF3 통합 1채널 + CMB + Cal = 3채널 -> 6.5sigma/7.2sigma")
print(f"  v23.2:  DR / CF3 / Cal = 3채널 (방법론 분리) + CMB 미포함 -> 5.5sigma/7.8sigma")
print(f"  ※ 9.4sigma: 이전 채널 구성 오류로 폐기. v23.2에서 미사용.")
print(f"  ※ DR(Cosmicflows-2, 위너필터) vs CF3(Cosmicflows-3, TF/SBF): 방법론적 독립")

print(f"\n[p_cf3 출처]")
print(f"  p_cf3 = 3.31e-4 -> validation4_CF3_regression.py 실행 결과")
print(f"  (c=+378.2+/-1.1 km/s 회귀의 t-검정 p값, N=15,174)")
