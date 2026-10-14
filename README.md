# rtlscout-bench-data

Run records and kept artefacts of the [RTL Scout leaderboard](https://github.com/plex1/rtlscout-bench). This
repository holds data only: no code lives here and nothing here imports anything. The tooling, the benchmark
suite and the page template are in [`rtlscout-bench`](https://github.com/plex1/rtlscout-bench), which includes
this repository as its `data/` submodule.

Everything in `runs/` and `campaigns/` is written by `python -m rtlscout_bench.record`; the two curation files
(`models.json`, `leaderboards/*.json`) are edited by hand or through that same tool.

**Scope: open-weights models only.** While `policy.open_weights_only` is on in the bench repository's
`campaign_matrix.yaml`, the recorder refuses every run of a model that is not marked `"open_weights": true` in
`models.json`, and the checks below fail if one is found.

## Layout

```
runs/<task>/<entry>/<run-id>/
    record.json        everything a table or a plot needs about the run (schema below)
    best_design/       the workspace as it was when the best design was found, plus _best_meta.json
    summary.txt        the agent's own write-up at the end of the run
    chat_log.txt       the full transcript: system prompt, tool calls, evaluation reports
campaigns/<name>.json  one file per campaign: versions, measured starting points, the list of runs
leaderboards/<task>.json   one per task: which runs count, exclusions with reasons, notes, starting points
models.json            display name, vendor, open weights, list prices, per model
```

- `<task>` and `<entry>` are the ids from `campaign_matrix.yaml`: a task is a benchmark with its cost metric, an
  entry is a model under one fixed set of agent conditions and is one leaderboard row.
- `<run-id>` is the time stamp rtlscout gave the run directory, `YYYYMMDD_HHMMSS`, in the clock of the machine
  that ran it (UTC in the published container image).
- Not stored: the per-evaluation workspaces, netlists and flow logs of a run. They are large and disposable.
  Whatever a later table or plot may need from them is therefore copied into `record.json` at recording time.

## `record.json`

One JSON object per run, one top-level key per line, one evaluation per line.

### Provenance

| field | type | meaning |
|---|---|---|
| `schema` | int | record format version, currently `1` |
| `suite` | string | suite tag the run was measured under, e.g. `suite-v1`; `pre-v1` for seed runs |
| `origin` | string | `campaign` (launched by `rtlscout_bench.campaign`), `seed` (older run adopted as initial data) or `manual` (other run directory adopted by hand) |
| `status` | string | `completed`, or `ended-early` when rtlscout reported an error for the run |
| `rtlscout` | string | version of the rtlscout package that produced the run, captured when the campaign was launched; an install from git adds its commit (`0.2.0+git.f4152b6`); free text for seed runs |
| `spire_hdl` | string or null | version of spire-hdl installed at launch; null when not recorded |
| `bench_commit` | string or null | commit of `rtlscout-bench` in use: at launch for campaign runs, at recording time for adopted runs; ends in `-dirty` if the checkout had uncommitted changes |
| `campaign` | string | campaign name, i.e. `campaigns/<name>.json`; only for `origin: campaign` |

### Conditions

Two runs belong in the same leaderboard row only if all of these fields and `suite` are identical.

| field | type | meaning |
|---|---|---|
| `task` | string | task id |
| `benchmark` | string | benchmark directory in `rtlscout-bench/benchmarks/`, e.g. `tpu_verilog` |
| `language` | string | `verilog` or `spirehdl` |
| `cost_metric` | string | rtlscout `--cost-metric`, e.g. `area_delay_product` |
| `target_delay` | number or null | rtlscout `--target-delay` in ps when the task sets one; null means rtlscout's default applied |
| `technology` | string | process technology of the flow, e.g. `asap7` |
| `entry` | string | entry id |
| `model` | string | model spec as given to rtlscout, `<provider>:<model>` |
| `open_weights` | bool or null | value of `models.json` for this model when the run was recorded |
| `reasoning_effort` | string | `default` (the provider's default), `low`, `medium`, `high`, or `off` for thinking disabled (`RTLSCOUT_REASONING_EFFORT`) |
| `reasoning_effort_default` | string | present when `models.json` named the vendor's default effort for this model at recording time (e.g. `max`); what a `default` row ran, shown in brackets in the row label like an explicit setting |
| `reasoning_effort_source` | string | present only when the run directory itself does not record the effort: where it is known from |
| `agent_backend` | string | `react`, rtlscout's built-in loop |
| `max_steps` | int | step budget of the run |
| `flags` | list of strings | extra rtlscout flags: the task's scoring flags followed by the entry's agent flags |
| `env` | object | environment variables the task sets, e.g. evaluation timeouts |

### Result

| field | type | meaning |
|---|---|---|
| `started` | string | start of the run, `YYYY-MM-DDTHH:MM`, from the run id |
| `best` | object or null | the best passing design: `cost`, `area` (µm²), `delay` (ps), `eval` and `step` at which it appeared, `file`, `td` (the `target_delay` the evaluation used, the agent's own or the task default; pass it as `--target-delay` to reproduce the cost) in `best_design/`. null if no evaluation passed |
| `n_evals` | int | number of evaluations the agent ran |
| `n_failed` | int | evaluations that did not pass the testbench or the flow |
| `steps_used` | int | agent steps taken (at most `max_steps`) |
| `runtime_min` | number | wall-clock duration of the run in minutes |
| `tokens` | object | `inp`, `out`, `cache_w`, `cache_r`: input, output, cache-write and cache-read tokens of the whole run |
| `price_per_mtok` | [number, number] or null | list price in USD per million input and output tokens used for `usd` |
| `price_cache_per_mtok` | [number, number] | cache-write and cache-read price; present only when the run used a cache |
| `usd` | number or null | API cost of the run: tokens times prices. Recomputable from the two fields above |
| `llm_providers` | object, optional | how many LLM calls each upstream provider served, e.g. `{"Z.AI": 61}` (rtlscout ce38834 and later) |
| `finish_reasons` | object, optional | how the LLM calls ended, e.g. `{"tool_calls": 59, "length": 2}`; `length` = the completion cap was hit |
| `error` | string | rtlscout's error text; present only for `ended-early` |

### `evals`: one object per evaluation, in order

| key | type | meaning |
|---|---|---|
| `e` | int | evaluation number, from 1 |
| `s` | int | agent step in which it ran |
| `ok` | bool | passed the testbench and the flow |
| `c` | number or null | cost; null for a failed evaluation |
| `area`, `delay`, `power` | number or null | cell area in µm², critical-path delay in ps, power as reported by the flow |
| `td` | number or null | target delay in ps the flow synthesised for |
| `ctx` | int | tokens in the model's context window at that point |
| `t_min` | number or null | minutes since the start of the run when the evaluation finished |
| `tok_in`, `tok_out` | int or null | input and output tokens used by the run up to this evaluation |
| `tok_cw`, `tok_cr` | int | cache-write and cache-read tokens up to this evaluation; present only when non-zero |

`t_min`, `tok_in` and `tok_out` come from fields rtlscout 0.2.0 writes for every evaluation (`elapsed_s`,
`cumulative_token_usage`). Runs made before that have no per-evaluation token counts (`null`), and their `t_min`
was reconstructed when they were recorded, from the modification time of each `eval_<n>/result.json` in the raw
run directory relative to the run's start.

## `leaderboards/<task>.json`

| field | meaning |
|---|---|
| `task` | task id; must equal the file name |
| `title`, `description`, `metric_label`, `unit`, `flow` | texts shown on the page for this task |
| `baselines` | per language: `cost`, `area`, `delay` of the unmodified starting point, with the `suite`, `rtlscout` version and `campaign` (or a `source` text) that measured it. This is the reference of the "vs. start" column |
| `runs` | per entry id: the run ids that count. One row of the table is exactly this list |
| `exclude` | runs that are in this repository but do not count: `entry`, `run`, `reason`. The reason is shown on the page |
| `notes` | per entry id: a note shown under the table |
| `footnotes` | general notes shown under the table |

`record` adds each new run to `runs`. A run that ended early is added to `exclude` instead, with the error as the
reason. When an entry gets its first run of a new suite, its runs of the previous suite move to `exclude` with the
reason "superseded": one row never mixes suites. To take a run out by hand:

```bash
python -m rtlscout_bench.record --exclude <task>/<entry>/<run-id> --reason "why"
```

## `models.json`

`models` maps a model spec (`<provider>:<model>`, exactly as in `campaign_matrix.yaml`) to:

| field | meaning |
|---|---|
| `name` | display name on the page |
| `vendor` | who made the model |
| `open_weights` | `true` only if the model's weights are openly available. Required for every model while the open-weights policy is on |
| `default_reasoning_effort` | optional: the reasoning effort the vendor applies when none is requested (e.g. `max`, `high`, `on`), as documented and checked at the time; copied into every record of the model and shown in brackets on default rows. Leave it out when the default is not known or not borne out by the runs |
| `price_per_mtok` | `input` and `output` list price in USD per million tokens; optional `cache_write` and `cache_read` (both default to the input price) |
| `price_note` | where the price applies, as free text |

## `campaigns/<name>.json`

Written by the campaign runner and copied here without local paths: `campaign` (name, suite, start and end time,
the arguments it was started with, rtlscout and spire-hdl versions, bench commit), `baselines` (the starting
point of each task and language as measured by `rtlscout.run_eval`: cost, area, delay, the command) and `runs`
(task, entry, model, repeat number, run id, status, exit code, best cost, the command, and whether the run was
recorded).

## Provenance rules

1. **A record is written once and not edited.** A number changes only by adding runs or by excluding runs, and
   an exclusion always carries a reason.
2. **Every record names what produced it:** suite tag, rtlscout version, spire-hdl version, bench commit. The
   conditions come from the campaign file written at launch, never from the state of a checkout at recording time.
3. **Runs that should not count stay in the repository** and are listed under `exclude`. Nothing is deleted to
   improve a row.
4. **Transcripts are published verbatim.** A run whose transcript or artefacts contain an API-key pattern
   (`sk-ant-`, `sk-or-`, `API_KEY=`) is refused by the recorder and can never be committed through it.
5. **The starting point is measured, not assumed:** the shipped starting design is evaluated once per campaign
   with the task's metric and flags. An agent run's first evaluation is not used as the reference.
6. **Seed data is marked as what it is.** The first twelve runs were made in July 2026, before the suite was
   fixed and before rtlscout was packaged. They carry `"suite": "pre-v1"`, `"origin": "seed"` and
   `"rtlscout": "pre-0.2.0 (unpackaged, 2026-07)"`, and the page labels their rows "measured before suite-v1".
   They were checked against the original 830-cycle vector file of the benchmark; suite-v1 ships a 6,782-cycle
   file whose first 830 cycles are identical. Their reasoning effort is not recorded in the run directories: the
   split of the two Kimi K3 entries (default effort, launched 2026-07-25; high effort, launched 2026-07-26)
   comes from the campaign notes, and those records say so in `reasoning_effort_source`.

## Checks

On every push and pull request, CI checks out `rtlscout-bench` and runs

```bash
python -m rtlscout_bench.record --check --data <this repository>
python -m rtlscout_bench.site --data <this repository>
```

which verifies: the schema of every record; that each record sits in the directory its own fields name; the
open-weights policy; that every run has its transcript, and its `best_design/` when it has a best design; that
no stored text contains an API-key pattern; that the leaderboard files only refer to runs that exist and give a
reason for each exclusion; and that the site builds from the selection.

## Licence

The records, summaries, transcripts and curation files in this repository are published under
[CC BY 4.0](LICENSE) (Felix Arnold). Attribution: "RTL Scout leaderboard data, https://github.com/plex1/rtlscout-bench-data".

## Third-party material

The design files under `runs/*/*/*/best_design/` derive from the benchmark starting points, which come from
LogikBench (MIT, Copyright (c) 2025 Zero ASIC Corporation). See `benchmarks/NOTICE` and `benchmarks/LICENSE` in
`rtlscout-bench`. The transcripts quote rtlscout's system prompt, which includes reference implementations with
their own licence headers.
