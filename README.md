# LA-axis: Leo-Aquarius Axis Anisotropy Analysis

**저자**: Kim, Gun-Sik (Yettumon)  
**소속**: 독립 연구자 | Seoul, KR  
**연락처**: nearcop3@gmail.com  
**라이선스**: CC BY 4.0  
**Zenodo DOI**: https://doi.org/10.5281/zenodo.20695683

---

## 개요

LA-axis(Leo-Aquarius Axis, RA=330°, Dec=0°)는 방사 분력이 소멸하는 관측자 기준 양방향 축으로, 물리 원리에서 선험적으로 정의된다. 본 리포지토리는 논문의 6개 검증 코드를 포함한다.

**핵심 결과:**
- 보수적 Fisher: p=1.77×10⁻⁸ → **5.5σ** (귀무=50%)
- 물리적 Fisher: p=3.69×10⁻¹⁵ → **7.8σ** (귀무=33.7%)

---

## 파일 구성

| 파일 | 논문 섹션 | 데이터 | 핵심 결과 |
|------|-----------|--------|-----------|
| `validation1_DR_angle.py` | 4.1절 | 불필요 | 7.5°, prob=0.0043 |
| `validation2_CMB_vector.py` | 4.2절 | 불필요 | 19.1°, prob=0.0275, 투영 94.4% |
| `validation3_calibrator_bias.py` | 4.3절 | 선택적 | 72.7%, p=4.11×10⁻⁵ |
| `validation4_CF3_regression.py` | 4.4절 | CF3 필수 | c=+378.2±1.1 km/s |
| `validation5_fisher_combined.py` | 6절 | 불필요 | 5.5σ / 7.8σ |
| `validation6_CF3_grid_scan.py` | 3.2절 | CF3 필수 | RA=332.6°±0.4° |

---

## 설치 및 실행

```bash
git clone https://github.com/yettumon/LA-axis.git
cd LA-axis
pip install -r requirements.txt

# 데이터 불필요한 코드 (즉시 실행 가능)
python validation1_DR_angle.py
python validation2_CMB_vector.py
python validation3_calibrator_bias.py
python validation5_fisher_combined.py

# CF3 데이터 필요
python validation4_CF3_regression.py
python validation6_CF3_grid_scan.py
```

---

## CF3 데이터

파일명: `CF3_all_distance.csv`  
출처: Cosmicflows-3 (Tully et al. 2016)  
공개 링크: https://edd.ifa.hawaii.edu/

**파싱 주의사항:**
- `skiprows=5, header=None`
- 좌표: `Glon/Glat` 사용 (`RAJ/DeJ`는 hms/dms 형식 → 신호 소멸)
- z 계산: `z = Vmod / c` (Vhel 아님)

---

## 관련 논문

| 논문 | DOI |
|------|-----|
| TRT (Tilted Ruler Theory) | https://doi.org/10.5281/zenodo.20572072 |
| LA-axis (본 논문) | https://doi.org/10.5281/zenodo.20695683 |
| IFM | https://doi.org/10.5281/zenodo.20038132 |

---

## 인용

```
Kim, G.S. (Yettumon), 2026. The Leo-Aquarius Axis (LA-axis):
A Priori Derivation of the Observer-Centric Principal Axis
Where Radial Force Vanishes in Cosmic Expansion.
Zenodo. https://doi.org/10.5281/zenodo.20695683
```

---

## 라이선스

Creative Commons Attribution 4.0 International (CC BY 4.0)  
https://creativecommons.org/licenses/by/4.0/
