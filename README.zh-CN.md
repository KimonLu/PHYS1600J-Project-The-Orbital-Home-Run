<p align="center">
  <img src="assets/background.png" alt="在月球上打棒球的宇航员" width="52%">
</p>

<h1 align="center">月球轨道全垒打</h1>
<p align="center"><strong>PHYS1600J 课程项目</strong></p>

<p align="center">
  <a href="README.md">English</a> | <strong>中文</strong>
</p>

如果在月球上击出一颗足够快的棒球，它能否进入低月球轨道、绕月飞行，并从击球手身后返回？本仓库通过解析轨道力学、考虑真实地形的高保真数值模型、正式项目报告和双语交互式轨迹求解器，给出了我们对 PHYS1600J Problem B 的回答。

- **求解器网站：** [Orbital Home Run Solver](https://kimonlu.github.io/PHYS1600J-Project-The-Orbital-Home-Run/)
- **项目报告：** [Physics1600J_Project.pdf](Project/docs/Physics1600J_Project.pdf)
- **代码仓库：** [KimonLu/PHYS1600J-Project-The-Orbital-Home-Run](https://github.com/KimonLu/PHYS1600J-Project-The-Orbital-Home-Run)

## 项目信息

- **课程：** PHYS1600J — Honors Physics
- **题目：** Problem B — The Orbital Home Run
- **报告全名：** *The Orbital Home Run — From a Measure-Zero Ideal Orbit to a Terrain-Aware Lunar Return*
- **小组：** Team 09
- **小组成员：** Kemeng Lu、Jiahao Yin、Yuxin Wu
- **指导教师：** Zijie Qu 教授
- **学院：** 上海交通大学 Global College
- **主要成果：** 科学报告、可复现的 Python 参考模型、版本化衍生数据与图像，以及双语浏览器求解器
- **参考模型：** LOLA LDEM64 月球地形和 GRAIL GRGM1200B 月球重力场
- **实现技术：** Python、TypeScript、Three.js、Vite、LaTeX
- **许可证：** [MIT](LICENSE)

## 项目背景

在地球上，全垒打主要是一段受大气影响的短程抛体运动；但在几乎没有大气的月球上，一颗速度足够高的棒球必须被视为绕月运行的质点。因此，计算月球“第一宇宙速度”只是问题的起点。

一次能够“绕月后击中自己”的月球全垒打必须同时满足四个彼此独立的条件：

1. 轨迹始终受月球引力束缚；
2. 棒球完成至少一周飞行，且途中不与月球表面或真实地形相交；
3. 棒球返回的是随月球自转而移动的击球点，而不仅是惯性系中的初始空间位置；
4. 在加入真实重力场、地形、数值误差与不确定性后，结论仍具有明确意义。

项目由精确的球对称二体模型逐步发展到三维、地形感知的数值传播模型。最终参考计算使用完整的 LOLA LDEM64 栅格地形、经收敛性检验后截断至 600 阶和 600 次的 GRGM1200B 重力场、匀速月球自转，以及短时间弧上固定几何的地球和太阳差分潮汐。

## 项目结果

- **理想速度尺度：** 在月球平均半径处，圆轨道速度为 **1.680 km/s**，逃逸速度为 **2.376 km/s**。
- **月面发射条件：** 对光滑、无自转的球形月球，从参考月面出发的有界轨迹只有在严格水平发射且 $1 \leq v/v_c < \sqrt{2}$ 时才不会立即穿入月面。因此，可行的“速度—角度”集合为零测集，且在结构上极不稳定。
- **月球自转影响：** 自转不会直接禁止返回，但会改变相遇时间。赤道顺行圆轨道将在会合周期 $2\pi/(n-\Omega)$ 后再次遇到移动的球场。
- **具体的地形感知返回解：** 600 阶边值解从 $5.4296875^\circ\text{N},\,201.3671875^\circ\text{E}$ 出发；球心位于**参考球面上方 30.000 km**，即**局部 LDEM64 地形上方 19.243 km**。
- **求得的初始状态：** 月面相对发射速度为 **1661.4559 m/s**，仰角为 **0.531611°**，从正北顺时针计的方位角为 **88.383063°**；计划返回时间为 **111.439 min**。
- **地形净空与数值残差：** 最小 LDEM64 地形净空为 **10.279 km**，定时边值问题残差为 **4.690 mm**。后者只是确定性模型内部的优化器残差，并不代表毫米级的现实物理精度。
- **物理意义：** 返回速度对应约 **200 kJ** 的动能。该结果排除了人类徒手击出这种轨迹的可能性，也没有证明相应发射或捕获装置在工程上可行。
- **网页模型交叉验证：** 与 Python 直接计算的 600 阶轨迹相比，浏览器重力图集在验证轨道上的位置差异为：15 km 高度处 **6.06 m**，30 km 高度处 **4.04 m**。

本项目在定义清晰的模型内给出了可复现的数值初始状态，但不宣称完成了导航级或现实物理认证。真正的月球实验还需要带日期的 SPICE 几何与星历链、局部地形测绘、重力场与发射状态协方差、太阳辐射压，以及面向完整不确定性集合的鲁棒约束修正。

## 求解器网站介绍

[双语静态求解器](https://kimonlu.github.io/PHYS1600J-Project-The-Orbital-Home-Run/)完全在浏览器内运行。用户可以输入月固坐标系中的发射地点、离地高度、月面相对速度、仰角、方位角、传播时长、返回容差、假定总位置误差上界，以及是否启用固定几何的地球与太阳潮汐项。

求解器提供：

- 中英文界面无刷新切换，且切换时保留当前计算结果；
- 返回、不确定返回、地形碰撞、逃逸和规定时限内不返回五类已验证示例；
- 惯性系三维轨迹、同步月球自转，以及 1×–1000× 物理时间播放；
- 等距圆柱投影下的月固地面轨迹；
- LDEM64 地形净空、首次碰撞坐标，以及离开后的轨迹与移动发射点之间的最近距离；
- 模型适用域提示和 RK4 步长加倍数值诊断。

网页按需加载无损 LDEM64 地形瓦片；重力模型由解析计算的中心项和基于 600 阶 GRGM1200B 生成的分高度层加速度修正图集组成。旋转坐标系动力学包含科里奥利项和离心项，并使用四阶 Runge–Kutta 步长加倍方法积分。该网站用于课程项目中的可复现实验与探索，不是导航、着陆或安全认证工具。

## 项目文件树

### `Project/` — 科学模型与报告

```text
Project/
├── README.md
├── docs/
│   ├── main.tex                       # 报告 LaTeX 源文件
│   ├── references.bib                 # 参考文献
│   ├── setting.cls                    # 报告格式
│   ├── project_statement.pdf          # 原始项目题目
│   └── Physics1600J_Project.pdf       # 最终报告
├── scripts/
│   ├── 01_ideal_models.py             # 解析二体轨道族
│   ├── 02_sensitivity.py              # 角度、高度与蒙特卡洛分析
│   ├── 03_rotation_and_resonance.py   # 月球自转与返回共振
│   ├── 04_realistic_perturbations.py  # 低阶重力与第三体摄动层级
│   ├── 05_terrain_envelope.py         # 保守地形包络
│   ├── 06_validation.py               # 中心引力模型收敛性验证
│   ├── 07_gravity_convergence.py      # GRGM 截断阶数选择
│   ├── 08_prepare_web_data.py         # 生成网页地形与重力数据
│   ├── 09_validate_web_solver.py      # 网页模型与直接轨迹交叉验证
│   ├── 10_high_fidelity_case.py       # 600 阶边值问题求解
│   ├── 11_case_sensitivity.py         # 雅可比收敛与敏感性
│   ├── 12_validate_web_terrain.py     # 无损地形瓦片验证
│   ├── 13_terrain_visualizations.py   # 地形与轨迹可视化
│   ├── 14_surface_feasibility.py      # 有界大圆地形可行性分析
│   ├── 15_high_fidelity_validation.py # 积分器与重力阶数验证
│   ├── 16_height_continuation.py      # 有界发射高度延拓
│   ├── general_solver.py              # 任意输入的通用求解器
│   ├── lunar_gravity.py               # 球谐月球重力模型
│   ├── lunar_terrain.py               # LDEM64 读取与插值
│   ├── orbital_home_run.py            # 共用轨道动力学
│   ├── plotting.py                    # 共用绘图工具
│   ├── generate_all.py                # 复现流程入口
│   ├── validate_results.py            # 已提交结果验证
│   ├── download_science_data.py       # NASA 数据下载与校验
│   └── requirements.txt
├── data/
│   ├── input/                         # 常数与外部数据清单
│   ├── external/                      # 下载的 NASA 数据，不纳入版本控制
│   └── output/                        # 衍生 CSV/JSON/NPZ 结果
└── figures/                           # 用于报告的图像
```

### `Web/` — 双语浏览器求解器

```text
Web/
├── index.html
├── package.json
├── package-lock.json
├── tsconfig.json
├── README.md
├── src/
│   ├── main.ts                        # 用户界面与结果展示
│   ├── solver.worker.ts               # 轨迹积分 Web Worker
│   ├── moon-view.ts                   # 交互式三维月球视图
│   ├── ground-track-view.ts           # 月固地面轨迹视图
│   ├── data-store.ts                  # 科学数据流式加载
│   ├── examples.ts                    # 已验证示例轨迹
│   ├── i18n.ts                        # 中英文文本
│   ├── types.ts                       # 共用 TypeScript 类型
│   └── style.css
└── public/
    ├── assets/                        # 月球视觉资源
    └── data/
        ├── terrain/                   # 无损 LDEM64 瓦片与清单
        └── gravity/                   # 600 阶重力图集瓦片与清单
```

## 如何复现

### 1. 克隆仓库并获取 Git LFS 数据

建议安装 Git、[Git LFS](https://git-lfs.com/)、Python 3.10 或更高版本，以及 Node.js 22。只有重新编译报告时才需要 TeX Live 2025。

```bash
git clone https://github.com/KimonLu/PHYS1600J-Project-The-Orbital-Home-Run.git
cd PHYS1600J-Project-The-Orbital-Home-Run
git lfs install
git lfs pull
```

### 2. 复现 `Project/` 中的计算

安装 Python 依赖：

```bash
cd Project
python -m pip install -r scripts/requirements.txt
```

轻量级解析模型、敏感性分析、自转、摄动和中心引力验证不需要下载大型外部数据：

```bash
python scripts/generate_all.py --profile core
```

如需运行高保真计算，先下载权威的 LOLA LDEM64 和 GRAIL GRGM1200B 数据。下载器会验证预期文件大小与 SHA-256：

```bash
python scripts/download_science_data.py
python scripts/generate_all.py --profile science
```

运行全部计算、报告图像和网页数据生成流程：

```bash
python scripts/generate_all.py --profile full
```

高阶重力场与全球地形步骤需要较长计算时间。如只需验证仓库中已提交的机器可读结果而不重新生成：

```bash
python scripts/generate_all.py --profile validate
```

标称返回解保存在 [`data/output/high_fidelity_case_summary.json`](Project/data/output/high_fidelity_case_summary.json)，采样轨迹保存在 [`data/output/high_fidelity_case_trajectory.csv`](Project/data/output/high_fidelity_case_trajectory.csv)。

### 3. 复现 `Project/` 中的可视化与报告

`core` 和 `science` 两个 profile 会分别生成对应图像，并写入 `Project/figures/`。下载科学数据后，也可以单独重新生成地形可视化：

```bash
python scripts/13_terrain_visualizations.py
```

重新编译报告：

```bash
cd docs
latexmk -pdf main.tex
```

### 4. 在本地运行 `Web/` 求解器

仅运行浏览器求解器时，仓库中由 Git LFS 管理的数据已经足够。从仓库根目录执行：

```bash
cd Web
npm ci
npm run dev
```

Vite 会在终端中输出本地访问地址。复现生产构建并在本地预览：

```bash
npm run build
npm run preview
```

如果希望从 NASA 权威源重新生成网页地形与重力数据，而不是直接使用仓库内的 LFS 对象：

```bash
cd ../Project
python scripts/download_science_data.py
python scripts/generate_all.py --profile web
```
