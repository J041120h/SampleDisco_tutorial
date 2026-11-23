#!/usr/bin/env python
"""Generate the SampleDisco config key reference (docs/tutorials/config_reference.md).

Usage:
    python gen_config_reference.py [--out config_reference.md] [--version-label X.Y.Z]

Everything in the reference is read from the installed (or PYTHONPATH) package:

* the keys, their order and their defaults come from
  ``inspect.signature(sampledisco.wrapper.wrapper.wrapper)``;
* the demo values and inline comments come from the packaged template
  ``sampledisco/config/config_demo.yaml`` (the file ``sampledisco --init-config`` writes);
* "consumed by" comes from a static (AST) walk of ``wrapper.py`` and of the per-modality
  wrappers it calls (``rna_wrapper``, ``atac_wrapper``, ``multiomics_wrapper``,
  ``downstream_analysis``). For each key the script follows the value one level into the
  sub-wrapper and lists every pipeline function/parameter it is passed to, plus the stage
  switches it sits under. A key that is never read, or is read only under a switch the
  caller hard-codes to ``False``, is flagged **dead**.

The only hand-written content is the GOTCHAS table below (short notes for keys whose
behaviour is surprising). Re-run the script after changing the package; do not edit the
generated page by hand.
"""
import argparse
import ast
import importlib
import inspect
import re
import sys
import textwrap
from collections import OrderedDict, defaultdict
from importlib.resources import files

SUB_WRAPPERS = ("rna_wrapper", "atac_wrapper", "multiomics_wrapper", "downstream_analysis")
GROUPS = (("Global", ""), ("RNA", "rna_"), ("ATAC", "atac_"), ("Multi-omics", "multiomics_"))
# Calls whose result is the (transformed) input value: data flows through them.
PASS_THROUGH = {"_coerce_sample_level_batch_col_list", "_first_batch_col_for_raisin",
                "_remap_units_to_modality", "resolve_embedding_keys", "_resolve_embedding_keys",
                "join", "get"}
# Wrapper helpers shown under the public function they forward to.
DISPLAY_ALIAS = {"_sample_adata_or_none": "build_sample_adata"}
# Book-keeping arguments that carry many keys at once (not a real consumer).
IGNORED_KWARGS = {"status_flags", "status_file_path"}
# More consumers than this: list function names only.
MAX_CONSUMERS = 8

# ---------------------------------------------------------------------------------------
# Hand-written notes. Each entry: (regex on the key, note). "{p}" expands to the key's
# modality prefix (rna_ / atac_ / multiomics_). Keep them short and verifiable.
# ---------------------------------------------------------------------------------------
GOTCHAS = [
    (r"^initialization$",
     "Default `true`. `true` deletes `{output_dir}/{rna,atac,multiomics}/result/` (a folder no "
     "stage writes) and rewrites `sys_log/main_process_status.json`; `false` reloads that JSON "
     "and prints it. Neither setting skips a stage: the stage flags alone decide what runs."),
    (r"^save_intermediate$",
     "Read only by the multi-omics pipeline (re-saves `preprocess/adata_sample.h5ad` after the "
     "Harmony pass and after the sample embedding). RNA and ATAC always write their files."),
    (r"^large_data_need_extra_memory$",
     "GPU only: switches the RAPIDS memory pool to managed (oversubscribable) memory. No effect on CPU."),
    (r"^continue_on_error$",
     "`false`: the first failing modality aborts the run. `true`: the error is stored as "
     "`results['<modality>_error']`, later modalities still run, and the CLI still exits 1."),
    (r"^random_state$",
     "Seeds the RNA/ATAC sample embedding and RAISIN only. Preprocessing and cell typing use a "
     "fixed seed 42; multi-omics uses `multiomics_random_state`."),
    (r"^use_gpu$",
     "Off unless set: runs Harmony, neighbors, Leiden and UMAP on the GPU through RAPIDS (the "
     "embedding's k-means stays on scikit-learn). Falls back to CPU with a warning if RAPIDS "
     "cannot initialise; `gpu_available` in `sys_log/main_process_status.json` records the "
     "backend. Results differ from a CPU run; see [GPU acceleration](gpu.md)."),
    (r"^run_rna_pipeline$", "Default `true`: set it to `false` explicitly for ATAC-only or multi-omics-only runs."),
    (r"^(rna|atac)_count_data_path$", "Raw counts in `.X`. Required whenever this pipeline runs, even when resuming."),
    (r"^(rna|atac)_adata_cell_path$",
     "Resume input: read when `{p}preprocessing: false` (default "
     "`{{output dir}}/preprocess/adata_preprocessed.h5ad`). With `{p}derive_sample_embedding: false` "
     "the embedding is taken from this file's `uns['X_DR_sample']`."),
    (r"^(rna|atac)_sample_meta_path$|^multiomics_(rna|atac)_sample_meta_file$",
     "Optional sample-level CSV, one row per sample, joined on the sample column. Its columns "
     "replace same-named `.obs` columns."),
    (r"^(rna|atac)_cell_meta_path$", "Optional cell-level CSV with a `barcode` column matching `obs_names`."),
    (r"^(rna|atac|multiomics)_sample_col$",
     "A column literally named `sample` is safest. With a sample-metadata CSV, this column is "
     "kept and a `sample` copy is added."),
    (r"^(rna|atac)_sample_level_batch_col$",
     "A string or a list. Two or more columns give multi-covariate sample-level Harmony. Harmony "
     "runs only with at least 8 samples and more than one batch value."),
    (r"^multiomics_batch_col$",
     "A single column name. It is scGLUE's `batch_key`, the `harmonize_xglue` batch, and the "
     "sample-level Harmony batch."),
    (r"^multiomics_celltype_col$",
     "Joint cell typing writes its labels here. With `multiomics_cell_type_cluster: false` the "
     "labels must already be in this column."),
    (r"^(rna|atac)_celltype_col$",
     "Cell typing always writes its labels to `obs['cell_type']`, and once it has run the "
     "embedding and downstream steps use that column. This key is the label source with "
     "`{p}existing_cell_types: true`, and the column the embedding reads when "
     "`{p}cell_type_cluster: false`."),
    (r"^(rna|atac)_existing_cell_types$",
     "Use labels already in `obs[{p}celltype_col]`; they are copied to `obs['cell_type']` (and "
     "merged down to `{p}n_target_cell_clusters` if that is smaller, with a warning)."),
    (r"^rna_min_cells$",
     "Used twice: genes detected in fewer cells are dropped, and samples with fewer cells are dropped."),
    (r"^(rna|atac|multiomics)_sample_embedding_pca_components$",
     "Capped at n_units - 1 (8 demo samples give 7 PCs)."),
    (r"^(rna|atac|multiomics)_sample_embedding_rmd_weight$",
     "Validated by the CLI: `\"equal\"` or a positive number."),
    (r"^(rna|atac|multiomics)_sample_embedding_batch_method$",
     "`\"harmony\"`, `\"linear\"` or `\"none\"`; any other value is rejected by the CLI and by "
     "`compute_sample_embedding`."),
    (r"^(rna|atac|multiomics)_sample_embedding_(use_clr|use_rmd|block_weights|rmd_weight)$",
     "Ignored when `{p}autotune_enable: true`."),
    (r"^(rna|atac|multiomics)_save_cell_adata_after_embedding$",
     "Optional: may be left out of a config (default `true`)."),
    (r"^(rna|atac|multiomics)_autotune_enable$",
     "Replaces `compute_sample_embedding` with `run_autotune` (CPU only). Without a batch column "
     "and without `{p}autotune_grouping_col` no search runs and alpha stays \"equal\"."),
    (r"^(rna|atac|multiomics)_grouping_columns$",
     "**Change for your data.** `null` becomes `[\"sev.level\"]`. Every column must exist in the "
     "sample metadata, or the visualization step raises `ValueError: Missing required columns`."),
    (r"^(rna|atac|multiomics)_trajectory_col$",
     "**Change for your data.** Default `sev.level`. Coerced to numeric for CCA: text labels "
     "become NaN and those samples are dropped."),
    (r"^(rna|atac|multiomics)_trajectory_visualization_label$",
     "**Change for your data.** `null` becomes `[\"sev.level\"]` (TSCAN plot colouring)."),
    (r"^(rna|atac|multiomics)_n_cca_pcs$",
     "Capped at the number of embedding PCs. Keep it well below the sample count: with PCs close "
     "to n_samples - 1 CCA fits any labelling perfectly."),
    (r"^(rna|atac|multiomics)_cca_pvalue$",
     "Permutation test with the same number of PCs as the observed score. When the PCs approach "
     "n_samples - 1 every permuted score is ~1: the step warns and the p-value is meaningless."),
    (r"^(rna|atac|multiomics)_tscan_origin$", "An integer: the index of the endpoint cluster to start from."),
    (r"^(rna|atac|multiomics)_tscan_pseudotime_mode$", "`\"rank\"` or `\"distance\"`; anything else raises."),
    (r"^(rna|atac|multiomics)_spline_order$",
     "A gene is fitted only if it is non-zero in at least order + splines + 2 samples (after "
     "automatic reduction). With 8 samples use 2; order 3 fits nothing."),
    (r"^(rna|atac|multiomics)_cluster_differential_gene_group_col$",
     "When set, the proportion test and RAISIN compare these groups **instead of** the k-means "
     "sample clusters (the demo sets `sev.level`)."),
    (r"^(rna|atac|multiomics)_cluster_dge$",
     "RAISIN: memory-heavy (scales with cores and features). Off in the demo."),
    (r"^(rna|atac|multiomics)_proportion_test$",
     "Groups: `{p}cluster_differential_gene_group_col` if set, else the k-means clusters "
     "(needs `{p}sample_cluster: true`)."),
    (r"^(rna|atac|multiomics)_summary_sample_csv_path$",
     "Written by the distance step and deleted after the trajectory step: always for the default "
     "path, and for a path you supply only if the run created it."),
    (r"^multiomics_(rna|atac)_file$",
     "Raw counts for scGLUE training: needed only with `multiomics_integration: true`."),
    (r"^multiomics_integrated_h5ad_path$",
     "Input only: read when `{{multiomics output}}/preprocess/adata_sample.h5ad` does not exist "
     "yet; every update goes to that output file."),
    (r"^(rna|atac|multiomics)_sample_distance_methods$",
     "SciPy metrics on the embedding (`euclidean`, `cosine`, `correlation`, … 12 in all) and "
     "`EMD`, `chi_square`, `jensen_shannon` on the cell-type composition; unknown names are "
     "skipped with a warning."),
    (r"^multiomics_use_gpu$", "GPU runs only if the global `use_gpu` is also `true`."),
    (r"^multiomics_random_state$",
     "Seeds the `harmonize_xglue` pass only; the multi-omics sample embedding uses its built-in seed 42."),
    (r"^multiomics_target_modality$",
     "When set, the embedding plots also include modality-highlight plots "
     "`<modality>_<column>_{expression,proportion}.png`, coloured by `multiomics_color_col` or else "
     "the first grouping column. `null` skips them."),
    (r"^multiomics_color_col$",
     "Colour column of the modality-highlight plots (see `multiomics_target_modality`)."),
    (r"^multiomics_visualization_grouping_column$",
     "Without it the multi-omics embedding plots are uncoloured."),
    (r"^multiomics_cell_type_cluster$",
     "Labels RNA cells by Leiden and transfers them to ATAC. Modality values must be `RNA` / `ATAC`."),
]


# ---------------------------------------------------------------------------------------
# AST helpers
# ---------------------------------------------------------------------------------------
def call_name(func):
    """Display name of a call target: `a.b.c(...)` -> 'c', `f(...)` -> 'f'."""
    if isinstance(func, ast.Name):
        return func.id
    if isinstance(func, ast.Attribute):
        return func.attr
    return None


def value_names(node):
    """Names read by an expression. Does not descend into calls, except pass-through ones."""
    out = set()

    def rec(n):
        if isinstance(n, ast.Name):
            out.add(n.id)
            return
        if isinstance(n, ast.Call):
            if call_name(n.func) in PASS_THROUGH:
                for a in n.args:
                    rec(a)
                for k in n.keywords:
                    rec(k.value)
            return
        if isinstance(n, (ast.Lambda, ast.ListComp, ast.DictComp, ast.SetComp, ast.GeneratorExp)):
            return
        for c in ast.iter_child_nodes(n):
            rec(c)

    rec(node)
    return out


def target_names(t):
    if isinstance(t, ast.Name):
        return [t.id]
    if isinstance(t, (ast.Tuple, ast.List)):
        return [x for e in t.elts for x in target_names(e)]
    return []


def source_stub(mod, attr):
    """Stand-in for a function whose module cannot be imported: name, module and
    parameter names read from the source file (no import)."""
    import importlib.util
    import types
    try:
        origin = importlib.util.find_spec(mod).origin
        tree = ast.parse(open(origin, encoding="utf-8").read())
    except Exception:
        return None
    for n in tree.body:
        if isinstance(n, ast.FunctionDef) and n.name == attr:
            a = n.args
            stub = types.SimpleNamespace(__name__=attr, __module__=mod)
            stub._sig_names = [x.arg for x in a.posonlyargs + a.args + a.kwonlyargs]
            return stub
    return None


class FuncAnalysis:
    """Data flow of a function's parameters into calls and `if` tests."""

    def __init__(self, fdef, module):
        self.fdef = fdef
        self.module = module
        args = fdef.args
        self.params = [a.arg for a in args.posonlyargs + args.args + args.kwonlyargs]
        self.imports = {}       # local name -> object (resolved lazily)
        self.aliases = {}       # local name -> (cpu_name, gpu_name, test_names) for IfExp aliases
        self.simple_alias = {}  # local name -> name, for `f = g`
        self._collect_imports()
        self.derive = self._derive()
        self.calls = []         # (call_node, display, obj, conds)
        self.tests = []         # (if_node, conds)
        self.refcount = defaultdict(int)
        self._count_refs()
        self._walk(fdef.body, [])

    # -- imports and aliases ------------------------------------------------------------
    def _collect_imports(self):
        for node in ast.walk(self.fdef):
            if isinstance(node, ast.ImportFrom):
                mod = node.module or ""
                if node.level:
                    base = self.module.__name__.rsplit(".", node.level)[0]
                    mod = f"{base}.{mod}" if mod else base
                for a in node.names:
                    self.imports[a.asname or a.name] = (mod, a.name)
            elif isinstance(node, ast.Assign) and len(node.targets) == 1 \
                    and isinstance(node.targets[0], ast.Name):
                v = node.value
                name = node.targets[0].id
                if isinstance(v, ast.IfExp) and isinstance(v.body, ast.Name) \
                        and isinstance(v.orelse, ast.Name):
                    self.aliases.setdefault(name, (v.orelse.id, v.body.id, value_names(v.test)))
                elif isinstance(v, ast.Name) and name not in self.aliases:
                    self.simple_alias.setdefault(name, v.id)

    def resolve(self, name):
        """Return (display_name, object or None) for a called name."""
        if name in self.aliases:
            name = self.aliases[name][0]
        name = self.simple_alias.get(name, name)
        obj = None
        if name in self.imports:
            mod, attr = self.imports[name]
            try:
                obj = getattr(importlib.import_module(mod), attr)
            except Exception:
                # e.g. scGLUE not installed: read the signature from the source instead
                obj = source_stub(mod, attr)
        else:
            obj = getattr(self.module, name, None)
        disp = getattr(obj, "__name__", None) or name
        return DISPLAY_ALIAS.get(disp, disp), obj

    # -- data flow ------------------------------------------------------------------------
    def _derive(self):
        derive = {p: {p} for p in self.params}
        assigns = []
        for node in ast.walk(self.fdef):
            if isinstance(node, ast.Assign):
                src = value_names(node.value)
                for t in node.targets:
                    if (isinstance(t, (ast.Tuple, ast.List)) and isinstance(node.value, (ast.Tuple, ast.List))
                            and len(t.elts) == len(node.value.elts)):
                        # `a, b = x, y`: pair element-wise, not all targets with all sources
                        for te, ve in zip(t.elts, node.value.elts):
                            for n in target_names(te):
                                assigns.append((n, value_names(ve)))
                        continue
                    for n in target_names(t):
                        assigns.append((n, src))
            elif isinstance(node, (ast.AugAssign, ast.AnnAssign)) and node.value is not None:
                for n in target_names(node.target):
                    assigns.append((n, value_names(node.value)))
            elif isinstance(node, ast.For):
                for n in target_names(node.target):
                    assigns.append((n, value_names(node.iter)))
        changed = True
        while changed:
            changed = False
            for tgt, src in assigns:
                new = set()
                for s in src:
                    new |= derive.get(s, set())
                cur = derive.setdefault(tgt, set())
                if not new <= cur:
                    cur |= new
                    changed = True
        return derive

    def params_of(self, names):
        out = set()
        for n in names:
            out |= self.derive.get(n, set())
        return out

    def _count_refs(self):
        for stmt in self.fdef.body:
            for node in ast.walk(stmt):
                if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load) and node.id in self.params:
                    self.refcount[node.id] += 1

    def _walk(self, stmts, conds):
        for s in stmts:
            self._visit(s, conds)

    def _visit(self, node, conds):
        if isinstance(node, ast.If):
            self.tests.append((node, conds))
            self._visit_expr(node.test, conds)
            self._walk(node.body, conds + [(node.test, True)])
            self._walk(node.orelse, conds + [(node.test, False)])
            return
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            return
        for c in ast.iter_child_nodes(node):
            if isinstance(c, ast.stmt):
                self._visit(c, conds)
            elif isinstance(c, ast.excepthandler):
                self._walk(c.body, conds)
            else:
                self._visit_expr(c, conds)

    def _visit_expr(self, node, conds):
        for n in ast.walk(node):
            if isinstance(n, ast.Call):
                name = call_name(n.func)
                if name is None or name in PASS_THROUGH:
                    continue
                disp, obj = self.resolve(name)
                self.calls.append((n, name, disp, obj, conds))

    # -- queries --------------------------------------------------------------------------
    def kw_params(self, call):
        """{callee_param: set(params)} for a call node, positional args named when possible."""
        out = {}
        name = call_name(call.func)
        _, obj = self.resolve(name)
        sig_names = []
        if getattr(obj, "_sig_names", None) is not None:
            sig_names = obj._sig_names
        elif obj is not None:
            try:
                sig_names = [p for p in inspect.signature(obj).parameters]
            except (TypeError, ValueError):
                sig_names = []
        for i, a in enumerate(call.args):
            ps = self.params_of(value_names(a))
            if ps:
                key = sig_names[i] if i < len(sig_names) else f"arg{i}"
                out.setdefault(key, set()).update(ps)
        for k in call.keywords:
            if k.arg is None:
                continue
            ps = self.params_of(value_names(k.value))
            if ps:
                out.setdefault(k.arg, set()).update(ps)
        return out

    def gated_calls(self, if_node, polarity=True):
        names = []
        for s in (if_node.body if polarity else if_node.orelse):
            for n in ast.walk(s):
                if isinstance(n, ast.Call):
                    nm = call_name(n.func)
                    if nm is None:
                        continue
                    disp, obj = self.resolve(nm)
                    mod = getattr(obj, "__module__", "") or ""
                    if ((mod.startswith("sampledisco") or nm in self.aliases)
                            and not disp.startswith("_") and disp not in names):
                        names.append(disp)
        return names


# ---------------------------------------------------------------------------------------
# Condition translation (callee parameter names -> config keys)
# ---------------------------------------------------------------------------------------
class Sub(ast.NodeTransformer):
    def __init__(self, mapping):
        self.mapping = mapping
        self.ok = True

    def visit_Name(self, node):
        if node.id in self.mapping:
            v = self.mapping[node.id]
            if isinstance(v, ast.AST):
                return ast.copy_location(v, node)
            return ast.copy_location(ast.Name(id=v, ctx=node.ctx), node)
        self.ok = False
        return node


def translate(test, polarity, mapping):
    """Return (text, constant_value_or_None) of a condition in config-key terms, or None."""
    sub = Sub(mapping)
    t = sub.visit(ast.parse(ast.unparse(test), mode="eval").body)
    if not sub.ok:
        return None
    const = None
    if isinstance(t, ast.Constant):
        const = bool(t.value) if polarity else (not bool(t.value))
    elif isinstance(t, ast.UnaryOp) and isinstance(t.op, ast.Not) and isinstance(t.operand, ast.Constant):
        const = (not bool(t.operand.value)) if polarity else bool(t.operand.value)
    text = ast.unparse(t)
    if not polarity:
        text = f"not {text}" if re.fullmatch(r"[\w.]+", text) else f"not ({text})"
    return text, const


def off_reasons(callee, conds, sconds):
    """Describe the hard-coded switches that make a use unreachable."""
    out = []
    for (test, _), c in zip(conds, sconds):
        if c is not None and c[1] is False:
            out.append(f"`{ast.unparse(test)}`, which `wrapper()` hard-codes to `False` "
                       f"in this pipeline's `{callee}()` call")
    return out


# ---------------------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------------------
def analyse(wmod):
    tree = ast.parse(inspect.getsource(wmod))
    fdefs = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}
    top = FuncAnalysis(fdefs["wrapper"], wmod)
    subs = {"downstream_analysis": FuncAnalysis(fdefs["downstream_analysis"], wmod)}
    for name in ("rna_wrapper", "atac_wrapper", "multiomics_wrapper"):
        m = importlib.import_module(f"sampledisco.wrapper.{name}")
        st = ast.parse(inspect.getsource(m))
        fd = next(n for n in st.body if isinstance(n, ast.FunctionDef) and n.name == name)
        subs[name] = FuncAnalysis(fd, m)

    info = {p: {"uses": [], "dead": None, "conds": []} for p in top.params}

    def top_mapping(call):
        """Callee param -> config key (str) or Constant node, for one sub-wrapper call."""
        mp = {}
        for k in call.keywords:
            if k.arg is None:
                continue
            if isinstance(k.value, ast.Constant):
                mp[k.arg] = ast.Constant(k.value.value)
            else:
                ps = top.params_of(value_names(k.value))
                if len(ps) == 1:
                    mp[k.arg] = next(iter(ps))
        return mp

    def top_cond(conds):
        out = []
        for test, pol in conds:
            names = value_names(test)
            if any(n.startswith("run_") and n.endswith("_pipeline") for n in names):
                continue
            mp = {}
            for n in names:
                ps = top.params_of({n})
                if len(ps) == 1:
                    mp[n] = next(iter(ps))
            tr = translate(test, pol, mp)
            out.append(tr)
        return out

    for call, name, disp, obj, conds in top.calls:
        kp = top.kw_params(call)
        if not kp:
            continue
        tconds = top_cond(conds)
        if name in SUB_WRAPPERS:
            sub = subs[name]
            mp = top_mapping(call)
            for cparam, keys in kp.items():
                if cparam in IGNORED_KWARGS:
                    continue
                for key in keys:
                    if cparam not in sub.params:
                        continue
                    info[key].setdefault("via", []).append((name, cparam, sub, mp, tconds))
        else:
            mod = getattr(obj, "__module__", "") or ""
            if not (mod.startswith("sampledisco") or name == "reinitialize") or disp.startswith("_"):
                continue
            for cparam, keys in kp.items():
                if cparam in IGNORED_KWARGS:
                    continue
                for key in keys:
                    info[key]["uses"].append((disp, cparam, tconds, False))

    # follow each key one level into the sub-wrappers
    for key, rec in info.items():
        for name, cparam, sub, mp, tconds in rec.get("via", []):
            if sub.refcount[cparam] == 0:
                rec.setdefault("dead_reasons", []).append(
                    f"`{name}()` accepts `{cparam}` but never reads it")
                continue
            found = False
            for call, cname, disp, obj, conds in sub.calls:
                kp = sub.kw_params(call)
                for callee_param, ps in kp.items():
                    if cparam not in ps or callee_param in IGNORED_KWARGS:
                        continue
                    mod = getattr(obj, "__module__", "") or ""
                    if not mod.startswith("sampledisco") or disp.startswith("_"):
                        continue
                    sconds = [translate(t, p, mp) for t, p in conds]
                    unreachable = any(c is not None and c[1] is False for c in sconds)
                    if unreachable:
                        rec.setdefault("off_reasons", []).extend(off_reasons(name, conds, sconds))
                    rec["uses"].append((disp, callee_param, tconds + sconds, unreachable))
                    found = True
            for if_node, conds in sub.tests:
                if cparam not in sub.params_of(value_names(if_node.test)):
                    continue
                sconds = [translate(t, p, mp) for t, p in conds]
                unreachable = any(c is not None and c[1] is False for c in sconds)
                if unreachable:
                    rec.setdefault("off_reasons", []).extend(off_reasons(name, conds, sconds))
                gated = sub.gated_calls(if_node, True) or sub.gated_calls(if_node, False)
                rec["uses"].append((f"switch in {name}", ", ".join(gated[:4]), tconds + sconds,
                                    unreachable))
                found = True
            # GPU/CPU function selection via an IfExp alias (e.g. preprocess_gpu if use_gpu else preprocess)
            for alias, (cpu, gpu, test_names) in sub.aliases.items():
                if cparam in sub.params_of(test_names):
                    rec["uses"].append((f"switch in {name}", f"{gpu} instead of {cpu}", tconds, False))
                    found = True
            if not found:
                rec["uses"].append((name, cparam, tconds, False))

        # top-level `if` tests in wrapper() itself
        for if_node, conds in top.tests:
            if key in top.params_of(value_names(if_node.test)):
                if key.startswith("run_") and key.endswith("_pipeline"):
                    gated = [n for n in top.gated_calls(if_node) if n in SUB_WRAPPERS]
                    rec["uses"].append(("switch in wrapper", ", ".join(gated), [], False))
                else:
                    rec["uses"].append(("switch in wrapper", "", top_cond(conds), False))

        if top.refcount[key] == 0:
            rec.setdefault("dead_reasons", []).append("`wrapper()` never reads it")

        live = [u for u in rec["uses"] if not u[3]]
        if rec["uses"] and not live and not rec.get("dead_reasons"):
            why = sorted(set(rec.get("off_reasons", []))) or ["a switch hard-coded to `False`"]
            rec["dead_reasons"] = ["only used under " + " and ".join(why)]
        rec["dead"] = bool(rec.get("dead_reasons")) and not live
        # conditions common to every live use
        cond_sets = []
        informative = [u for u in live if not (u[0].startswith("switch in ") and not u[1])]
        for u in (informative or live):
            cond_sets.append({c[0] for c in u[2] if c is not None})
        rec["conds"] = sorted(set.intersection(*cond_sets)) if cond_sets else []
    return info


# ---------------------------------------------------------------------------------------
# Template parsing
# ---------------------------------------------------------------------------------------
def read_template():
    import yaml
    text = (files("sampledisco") / "config" / "config_demo.yaml").read_text()
    values = yaml.safe_load(text)
    comments = {}
    for line in text.splitlines():
        m = re.match(r"^([A-Za-z_]\w*):(.*)$", line)
        if not m:
            continue
        rest = m.group(2)
        # strip a trailing comment that is not inside quotes
        q, cpos = None, None
        for i, ch in enumerate(rest):
            if ch in "\"'":
                q = None if q == ch else (ch if q is None else q)
            elif ch == "#" and q is None and (i == 0 or rest[i - 1].isspace()):
                cpos = i
                break
        if cpos is not None:
            comments[m.group(1)] = rest[cpos + 1:].strip()
    return text, values, comments


def signature_sections(wmod):
    """Section label for each wrapper() parameter, from the comments in its signature."""
    src = inspect.getsource(wmod)
    tree = ast.parse(src)
    fdef = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "wrapper")
    lines = src.splitlines()
    args = fdef.args.args + fdef.args.kwonlyargs
    out, label, prev = {}, "Output", fdef.lineno
    for a in args:
        block = []
        for ln in range(prev, a.lineno - 1):
            s = lines[ln].strip()
            if s.startswith("#"):
                block.append(s.lstrip("#").strip())
            elif s and block and not s.startswith("#"):
                block = []
        cand = [b.strip("= ").strip() for b in block if b.strip("= ")]
        if any(c.endswith("PIPELINE PARAMETERS") for c in cand):
            label = "Input and output paths"
        elif cand:
            first = re.sub(r"\s*\(.*$", "", cand[0]).split(" -- ")[0]
            first = re.sub(r"^V2\s+", "", first).split(". ")[0].rstrip(".:")
            if (len(first) <= 48 and "=" not in first and ":" not in first
                    and not first.lower().startswith("todo")):
                label = first[0].upper() + first[1:]
        out[a.arg] = label
        prev = a.lineno
    return out


# ---------------------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------------------
def fmt(v):
    if v is inspect.Parameter.empty:
        return "*(required)*"
    if v is None:
        return "null"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, str):
        return f'"{v}"'
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(fmt(x) for x in v) + "]"
    if isinstance(v, float) and v >= 1000 and v == int(v):
        return f"{v:.1f}"
    return repr(v)


def norm(v):
    if isinstance(v, tuple):
        return [norm(x) for x in v]
    if isinstance(v, list):
        return [norm(x) for x in v]
    return v


def cell(s):
    return str(s).replace("|", "\\|").replace("\n", " ")


def gotchas_for(key):
    prefix = next((p for _, p in GROUPS if p and key.startswith(p)), "")
    out = []
    for pat, note in GOTCHAS:
        if re.search(pat, key):
            out.append(note.replace("{p}", prefix).replace("{{", "{").replace("}}", "}"))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default="-")
    ap.add_argument("--version-label", default=None)
    a = ap.parse_args()

    import sampledisco
    wmod = importlib.import_module("sampledisco.wrapper.wrapper")
    from sampledisco.cli import OPTIONAL_KEYS
    sig = inspect.signature(wmod.wrapper)
    params = list(sig.parameters.values())
    _, tvals, tcomments = read_template()
    sections = signature_sections(wmod)
    info = analyse(wmod)
    version = a.version_label or getattr(sampledisco, "__version__", "?")

    keys = [p.name for p in params]
    missing_in_template = [k for k in keys if k not in tvals and k not in OPTIONAL_KEYS]
    extra_in_template = [k for k in tvals if k not in sig.parameters]
    dead = [k for k in keys if info[k]["dead"]]

    L = []
    w = L.append
    w("# Config key reference")
    w("")
    w("<!-- GENERATED FILE: do not edit. Regenerate with docs/assets/gen_config_reference.py -->")
    w("")
    w(f"Every key that `sampledisco -m complex --config` accepts, for sampledisco **{version}**: "
      f"**{len(keys)} keys** "
      + "(" + ", ".join(f"{g} {sum(1 for k in keys if group_of(k) == g)}" for g, _ in GROUPS) + ").")
    w("")
    w("Generated from the package by [`gen_config_reference.py`](../assets/gen_config_reference.py); "
      "do not edit by hand. **Default** is the `wrapper()` default, which a key left out of a Python "
      "call takes; the CLI requires every key except the "
      f"{len(OPTIONAL_KEYS)} optional ones. **Demo** is the value in the template that "
      "`sampledisco --init-config` writes, in bold when it differs. **Consumed by** is the function "
      "and parameter the value reaches, or the stage it switches (`switch in X`). How the keys fit "
      "together is in the [Configuration guide](configuration.md).")
    w("")
    w(f'!!! warning "{len(dead)} keys are accepted but ignored"')
    w("    The CLI requires them, so leave them in the file, but changing them has no effect:")
    w("")
    for k in dead:
        w(f"    - `{k}`: {'; '.join(info[k]['dead_reasons'])}.")
    w("")
    w('??? info "Template vs. signature check"')
    w(f"    - Keys in `wrapper()` but missing from the template (excluding optional keys): "
      f"{', '.join(f'`{k}`' for k in missing_in_template) or 'none'}.")
    w(f"    - Keys in the template that `wrapper()` does not accept: "
      f"{', '.join(f'`{k}`' for k in extra_in_template) or 'none'}.")
    w(f"    - Optional keys (the CLI fills in the default): {', '.join(f'`{k}`' for k in sorted(OPTIONAL_KEYS))}.")
    w("")

    for gname, prefix in GROUPS:
        gkeys = [k for k in keys if group_of(k) == gname]
        w(f"## {gname} keys ({len(gkeys)})")
        w("")
        by_sec = OrderedDict()
        for k in gkeys:
            by_sec.setdefault(sections.get(k, "General"), []).append(k)
        for sec, ks in by_sec.items():
            w(f"### {gname}: {sec[0].lower() + sec[1:] if sec[:2] != sec[:2].upper() else sec}")
            w("")
            w("| Key | Default | Demo | Consumed by | Notes |")
            w("|---|---|---|---|---|")
            for k in ks:
                p = sig.parameters[k]
                d = p.default
                if k in tvals:
                    tv = tvals[k]
                    demo = fmt(tv)
                    if norm(tv) != norm(d):
                        demo = f"**{demo}**"
                else:
                    demo = "*(not in template)*"
                rec = info[k]
                live = [u for u in rec["uses"] if not u[3]]
                informative = [u for u in live if not (u[0].startswith("switch in ") and not u[1])]
                cons, funcs = [], []
                for disp, cp, _, _ in (informative or live):
                    if disp.startswith("switch in "):
                        txt = f"{disp}" + (f": {cp}" if cp else "")
                    else:
                        txt = f"`{disp}:{cp}`"
                        if f"`{disp}`" not in funcs:
                            funcs.append(f"`{disp}`")
                    if txt not in cons:
                        cons.append(txt)
                if len(cons) > MAX_CONSUMERS:
                    cons = [f"{len(funcs)} functions: " + ", ".join(funcs)]
                if rec["dead"]:
                    cons = ["*nothing*"]
                notes = []
                if rec["dead"]:
                    notes.append("**Dead key:** " + "; ".join(rec["dead_reasons"]) + ".")
                if rec["conds"] and not rec["dead"]:
                    notes.append("Only used when " + " and ".join(f"`{c}`" for c in rec["conds"]) + ".")
                notes += gotchas_for(k)
                if k in tcomments and tcomments[k]:
                    notes.append(f"Template: *{tcomments[k]}*")
                w("| " + " | ".join(cell(x) for x in (
                    f"`{k}`", f"`{fmt(d)}`" if d is not inspect.Parameter.empty else fmt(d),
                    (f"**`{demo[2:-2]}`**" if demo.startswith("**")
                     else demo if demo.startswith("*(") else f"`{demo}`"),
                    "<br>".join(cons) or "–",
                    " ".join(notes) or "")) + " |")
            w("")

    text = "\n".join(L) + "\n"
    if a.out == "-":
        sys.stdout.write(text)
    else:
        with open(a.out, "w", encoding="utf-8") as fh:
            fh.write(text)
        print(f"wrote {a.out}: {len(keys)} keys, {len(dead)} dead", file=sys.stderr)


def group_of(key):
    for g, p in GROUPS[1:]:
        if key.startswith(p):
            return g
    return "Global"


if __name__ == "__main__":
    main()
