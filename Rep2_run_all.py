"""
Rep2 - 一键复现 Rep1 的所有关键分析
依次运行:
  1. cellbin / bin20 / nucleus 三种分析 (不区分 IGHKL)
  2. 三者整合折线图
  3. cell segmentation 6 张随机 ROI 图
  4. ROI 中心区 BCR 叠加图 (DAPI + cellbin / nucleus)
  5. ROI 在全局的位置标注图
  6. BCR in/out 三口径判断图
  7. cellbin IGH×IGK 热图 (targetSeq + coords)
"""
import os, subprocess, sys

REP2_DIR = "/data/scratch/projects/punim1236/AAA_chen_here/AAA克隆_XCR/Rep2"
SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))

steps = [
    "Rep2_step1_three_QC.py",        # cellbin/bin20/nucleus 不区分 IGHKL
    "Rep2_step2_integrated_line.py", # 整合折线图
    "Rep2_step3_ROI_visuals.py",     # cell seg QC + BCR 叠加 + in/out
    "Rep2_step4_IGHK_heatmap.py",    # cellbin 的 IGH×IGK 热图
]

for s in steps:
    print(f"\n{'='*70}\nRunning {s}\n{'='*70}")
    rc = subprocess.run([sys.executable, os.path.join(SCRIPTS_DIR, s)])
    if rc.returncode != 0:
        print(f"FAILED: {s}")
        sys.exit(rc.returncode)

print("\n\nALL DONE. Outputs in:", REP2_DIR)