# DNA 限制性酶切模拟 API

基于 **Python + FastAPI + SQLite** 的限制性内切酶消化模拟服务。提交线性或环形 DNA
序列（仅 IUPAC 碱基）以及多种酶的识别模式与正/反链切口偏移，服务在**两个方向**
查找识别位点（允许重叠、含简并碱基、可跨越环形序列首尾），对同一物理切口去重后
按拓扑结构生成片段，返回坐标、长度、两端的平端/黏端类型及突出端序列。

## 快速启动

一条指令完成构建、依赖安装与数据库初始化：

```bash
docker compose up --build
```

依赖在镜像构建时安装（`pip install -r requirements.txt`），SQLite 数据库在容器
启动时自动建库建表，无需任何手工步骤。

也可以不使用 compose，直接用 Docker 构建运行：

```bash
docker build -t digest-api .
docker run -d --name digest-api -p 8000:8000 -v digest-data:/app/data digest-api
```

启动后验证服务可用：

```bash
curl http://localhost:8000/health   # {"status":"ok"}
```

停止服务：`docker compose down`（或 `docker stop digest-api`）。数据库保存在
`digest-data` 数据卷中，重启后数据仍在；`docker compose down -v` 会连同数据
一并删除。

## 配置

| 配置项 | 默认值 | 说明 |
| --- | --- | --- |
| 环境变量 `DATABASE_PATH` | `/app/data/digest.db`（容器内） | SQLite 数据库文件路径，目录不存在时自动创建 |
| 端口映射 | `8000:8000`（compose 中 `ports`） | 宿主机端口：容器端口 |
| 数据卷 | `digest-data` → `/app/data` | 持久化数据库文件 |

不使用 Docker 时本地运行：

```bash
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## API 地址

- 服务地址：**http://localhost:8000**
- 交互式文档（Swagger UI）：**http://localhost:8000/docs**
- 健康检查：`GET /health`

### 接口一览

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| `POST` | `/digests` | 提交序列与酶列表，执行酶切、保存并返回完整结果（201） |
| `GET` | `/digests` | 列出全部提交的摘要（按 id 升序） |
| `GET` | `/digests/{id}` | 按 id 查询保存的完整酶切结果（不存在返回 404） |
| `POST` | `/digests/{id}/variant-analysis` | 对已保存记录应用单碱基替换并重跑酶切，返回差异分析（不落库） |
| `GET` | `/health` | 健康检查 |

### 提交示例

```bash
curl -X POST http://localhost:8000/digests \
  -H 'Content-Type: application/json' \
  -d '{
    "sequence": "AAAAGAATTCTTTT",
    "topology": "linear",
    "enzymes": [
      {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5},
      {"name": "EcoRV", "pattern": "GATATC", "cut_top": 3, "cut_bottom": 3}
    ]
  }'
```

- `sequence`：仅允许 IUPAC 碱基 `ACGTRYSWKMBDHVN`（大小写均可，存储为大写）。
- `topology`：`linear`（线性）或 `circular`（环形）。
- `enzymes`：至少一种酶。
  - `pattern`：识别模式，同样仅允许 IUPAC 碱基；
  - `cut_top` / `cut_bottom`：正/反链切口偏移，以识别位点 5' 端为基准，
    合法范围为 **`0 ~ len(pattern)`**。

响应（节选）：

```json
{
  "id": 1,
  "topology": "linear",
  "sequence_length": 14,
  "cuts": [
    {
      "top_cut": 5,
      "bottom_cut": 9,
      "end_type": "5prime",
      "overhang": "AATT",
      "enzymes": ["EcoRI"],
      "sites": [{"enzyme": "EcoRI", "start": 4, "strand": "+", "matched": "GAATTC"}]
    }
  ],
  "fragments": [
    {"index": 0, "start": 0, "end": 5, "length": 5, "sequence": "AAAAG",
     "wraps_origin": false,
     "left_end": {"type": "none", "overhang": ""},
     "right_end": {"type": "5prime", "overhang": "AATT"}},
    {"index": 1, "start": 5, "end": 14, "length": 9, "sequence": "AATTCTTTT",
     "wraps_origin": false,
     "left_end": {"type": "5prime", "overhang": "AATT"},
     "right_end": {"type": "none", "overhang": ""}}
  ]
}
```

## 变异影响分析

`POST /digests/{id}/variant-analysis` 对已保存的记录应用一组 **0-based 单碱基
替换（SNV）**，在内存中对突变后的序列重跑酶切并与原结果对比。**分析只读取记录，
不写数据库**——原记录保持不变，也不会产生新记录。

```bash
curl -X POST http://localhost:8000/digests/1/variant-analysis \
  -H 'Content-Type: application/json' \
  -d '{"variants": [{"position": 9, "base": "T"}]}'
```

- `variants`：至少 1 个替换；`position` 为 0-based 替换位置，`base` 为单个
  IUPAC 碱基（大小写均可，存储为大写）。
- 校验（均返回 **422**，且不落库）：`position` 越界（`< 0` 或 `≥ 序列长度`）、
  `base` 非法（非单个 IUPAC 碱基）、同一 `position` 重复出现、`variants` 为空；
  记录不存在返回 **404**。

响应（节选，基于上文 `AAAAGAATTCTTTT` + EcoRI 的记录）：

```json
{
  "submission_id": 1,
  "topology": "linear",
  "sequence_length": 14,
  "variants_applied": [{"position": 9, "ref_base": "C", "base": "T"}],
  "mutated_sequence": "AAAAGAATTTTTTT",
  "added_cuts": [],
  "removed_cuts": [
    {"top_cut": 5, "bottom_cut": 9, "end_type": "5prime", "overhang": "AATT",
     "enzymes": ["EcoRI"],
     "sites": [{"enzyme": "EcoRI", "start": 4, "strand": "+", "matched": "GAATTC"}]}
  ],
  "changed_cuts": [],
  "fragment_lengths": {
    "baseline_lengths": [5, 9],
    "variant_lengths": [14],
    "added_lengths": [14],
    "removed_lengths": [5, 9]
  }
}
```

字段说明：

- **`variants_applied`**：实际应用的替换（含原碱基 `ref_base`），按位置升序。
- **`added_cuts` / `removed_cuts`**：变异后新增 / 消失的完整切口记录，
  结构同 `POST /digests` 的 `cuts`。
- **`changed_cuts`**：坐标不变但末端特征变化的切口，`changes` 列出发生变化的
  字段（`end_type` / `overhang`），`before` / `after` 为变化前后的末端。
  典型场景：替换落在黏端突出序列内部且未破坏识别位点（如模式含简并碱基）。
- **`fragment_lengths`**：两次酶切的片段长度序列（按片段顺序）及多重集差
  （`added_lengths` / `removed_lengths`，升序）。
- **坐标约定**：单碱基替换不改变序列长度，同一切口的坐标不会平移；识别位点的
  破坏/新建体现为一条 `removed` 加一条 `added`（即"坐标改变"）。
- **稳定排序**：切口列表按 `(top_cut, bottom_cut)` 排序，`variants_applied`
  按位置排序，长度差升序——相同输入的响应完全一致。

## 结果字段约定

- **坐标**：均为 0-based，作用于 top 链（提交序列方向）；切口位置 `c` 表示断裂
  位于第 `c-1` 与第 `c` 个碱基之间。片段区间为半开区间 `[start, end)`。
- **末端类型 `end_type` / `left_end.type` / `right_end.type`**：
  - `blunt`：平端（正/反链切口位置相同）；
  - `5prime`：5' 黏端；`3prime`：3' 黏端；
  - `none`：线性分子的天然末端，或未切割的完整分子。
- **`overhang`**：突出端序列，按 top 链 5'→3' 方向取自两切口之间的区段；
  平端与 `none` 时为空字符串。
- **环形序列**：切口位置对序列长度取模；`wraps_origin: true` 表示该片段跨越
  首尾原点，其 `end` 为取模后的坐标（可能小于 `start`）。单个切口开环时返回
  一条全长片段；无切口时保留整条环状分子。
- **结果顺序**：切口按 `(top_cut, bottom_cut)` 排序，片段按起点排序，酶名与
  位点列表均排序——**相同输入的结果顺序必然一致**。

## 行为规则

- **双向搜索**：同时在 top 链搜索 `pattern` 与 `reverse_complement(pattern)`
  （即酶结合于反链的情形）；回文模式只搜索一次。
- **重叠位点**：允许识别位点相互重叠，每个命中都会计算切口。
- **简并碱基**：模式与序列中的 IUPAC 简并碱基按"可能集合有交集"匹配
  （如模式中的 `A` 可匹配序列中的 `R`）。
- **跨原点**：环形序列在扩展序列上搜索，跨越首尾的识别位点与切口均可命中。
- **物理切口去重**：不同酶/不同位点若产生相同的 `(top_cut, bottom_cut)`，
  只保留一条切口记录，酶名合并进 `enzymes` 列表。
- **非法偏移整体拒绝**：任一酶的切口偏移超出 `[0, len(pattern)]`（含负数），
  或序列/模式含非 IUPAC 字符，整个请求以 **422** 拒绝，**不会有任何数据落库**。
- **变异分析不落库**：`POST /digests/{id}/variant-analysis` 只在内存中对突变
  序列重跑酶切并返回差异，原记录与数据库均不受影响。
- **无切口**：返回一条覆盖整条分子的片段，两端类型为 `none`。

## 项目结构

```
app/
├── main.py      # FastAPI 入口与路由
├── digest.py    # 酶切核心算法（搜索、去重、片段生成）
├── variant.py   # 变异影响分析（应用 SNV、切口/片段长度对比）
├── iupac.py     # IUPAC 编码表与反向互补
├── schemas.py   # 请求/响应模型与校验
└── db.py        # SQLite 持久化（启动时自动建表）
tests/           # pytest 测试（算法 + API + 变异分析）
Dockerfile
docker-compose.yml
requirements.txt
```

## 运行测试

```bash
pip install -r requirements-dev.txt
pytest
```
