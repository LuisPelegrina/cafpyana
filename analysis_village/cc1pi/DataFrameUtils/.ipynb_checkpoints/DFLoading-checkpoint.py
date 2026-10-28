import pandas as pd
import pyanalib.pandas_helpers as ph
import pyanalib.split_df_helpers as splh
import pyanalib.stat_helpers as sh
from analysis_village.cc1pi.DataFrameUtils import DFCleaning
from analysis_village.cc1pi.CutMasks import CutMasks
from analysis_village.cc1pi.Constants import CTE

import numpy as np

import gc
import numpy as np
import pandas as pd
from tqdm import tqdm
import sys
# Assumes the rest of DFLoading.py is available: splh, DFCleaning, _load_split, slielm_dfs


# ======================================================================
# shrink_df
# ======================================================================
# Columns that must never be touched (identifiers, POT sums, ...). A column is skipped
# if ANY element of its label matches (labels can be tuples for MultiIndex columns).
SHRINK_SKIP_COLS = {"run", "subrun", "evt", "event", "pot", "spill", "ntuple", "entry"}

_F32_MAX = 3.0e38
_F32_EXACT_INT = 2 ** 24      # float32 represents integers exactly only up to here


def _col_tokens(label):
    """Set of string tokens of a column label (works for tuple labels)."""
    if isinstance(label, tuple):
        return {str(t) for t in label}
    return {str(label)}

import gc
import pyarrow as pa
import pandas as pd

def zero_copy_concat_shift(dfs, level=0):
    """
    Concatenates DataFrames via PyArrow zero-copy tables, shifting the index level
    without allocating 2x memory during concatenation.
    """
    tables = []
    offset = 0
    idx_names = dfs[0].index.names if isinstance(dfs[0].index, pd.MultiIndex) else [dfs[0].index.name]
    target_idx = idx_names[level] if isinstance(level, int) else level

    # Process each DataFrame sequentially and immediately free the original Pandas object
    while dfs:
        df = dfs.pop(0)
        if df is None or df.empty:
            continue

        # Extract max value for shift offset
        max_val = int(df.index.get_level_values(target_idx).max())

        # Convert index to columns to allow Arrow zero-copy handling
        df_reset = df.reset_index()
        del df  # Free original reference immediately
        
        if offset > 0:
            df_reset[target_idx] += offset

        offset += max_val + 1

        # Convert to PyArrow Table
        tables.append(pa.Table.from_pandas(df_reset))
        del df_reset
        gc.collect()

    if not tables:
        return pd.DataFrame()

    # Zero-copy concatenation (allocates ~0 extra RAM; creates references to chunks)
    combined_table = pa.concat_tables(tables, promote_options="permissive")
    del tables
    gc.collect()

    final_df = combined_table.to_pandas(self_destruct=True, split_blocks=True)
    del combined_table; gc.collect()
    final_df.set_index(idx_names, inplace=True)
    return final_df

    
def _float_is_safe(s):
    """True if casting the float64 Series to float32 loses nothing important:
    no overflow, and no large integer-valued data (ids, counters) stored as float."""
    a = s.to_numpy()
    fin = a[np.isfinite(a)]
    if fin.size == 0:
        return True
    amax = np.abs(fin).max()
    if amax >= _F32_MAX:
        return False
    if amax >= _F32_EXACT_INT and np.all(fin == np.round(fin)):
        return False
    return True

def shrink_df(df, floats=True, ints=True, skip_cols=None, verbose=False, **_ignored):
    """In place, column by column:
       float64 -> float32  (only if safe: no overflow, no large integer-valued ids stored as float)
       int64   -> int32    (only if all values fit in int32)
    Everything else (objects, bools, categoricals, already-small dtypes) is untouched.
    The index is never modified. Do NOT use on the hdr DataFrame (POT sums need float64)."""
    if df is None or df.empty:
        return df

    skip = set(SHRINK_SKIP_COLS) | set(skip_cols or ())
    before = df.memory_usage(deep=True).sum() if verbose else 0
    i32 = np.iinfo(np.int32)

    for c in list(df.columns):
        if _col_tokens(c) & skip:
            continue
        s = df[c]
        if not isinstance(s, pd.Series):          # duplicated column label -> leave alone
            continue
        dt = s.dtype

        if floats and dt == np.float64:
            if _float_is_safe(s):
                df[c] = s.astype(np.float32)
        elif (ints and pd.api.types.is_integer_dtype(dt) and dt != np.bool_
              and np.dtype(dt).itemsize > 4 and len(s) > 0):
            if s.min() >= i32.min and s.max() <= i32.max:
                df[c] = s.astype(np.int32)

    if verbose:
        after = df.memory_usage(deep=True).sum()
        print(f"    shrink_df: {before / 1e9:.3f} GB -> {after / 1e9:.3f} GB")
    return df

import gc
import pandas as pd


def _apply_obvious_cosmic_filter(cc1pi_df):
    """Filters out obvious cosmic slices from cc1pi DataFrame using 6-level tuple key."""
    if cc1pi_df is None or cc1pi_df.empty:
        return cc1pi_df
    else:
        return cc1pi_df[cc1pi_df.slc.cut.obvious_cosmic]

# ======================================================================
# load_df
# ======================================================================
def load_df(file, keys2load, n_max_concat=100, filter_df=True,
            reprocess_df=True, reprocess_truth=True, shrink=True,
            cosmic_filter=True):

    print(f"keys in {file}")
    splh.print_keys(file)
    n_split = splh.get_n_split(file)
    print("n_split: %d" % n_split)
    n_load = min(n_split, n_max_concat)

    # ---- Pass 1: hdr only (small) -> global duplicate filtering ----
    hdr_clean = None
    if "hdr" in keys2load or filter_df:
        hdr_all = pd.concat(
            [_load_split(file, ["hdr"], i)["hdr"] for i in tqdm(range(n_load), desc="hdr")]
        )
        print(hdr_all["pot"].sum())
        if filter_df:
            print(f"duplication for {file}")
            DFCleaning.find_duplicate_run_evt_combinations(hdr_all)
            DFCleaning.plot_duplicate_run_subrun_evt_distribution(hdr_all, file)
            hdr_clean = DFCleaning.filter_unique_events(hdr_all)
            DFCleaning.find_duplicate_run_evt_combinations(hdr_clean)
        else:
            hdr_clean = hdr_all
        del hdr_all
        gc.collect()

    # ---- Pass 2: per split -> filter -> cosmic cut -> float32 -> keep ----
    other_keys = [k for k in keys2load if k != "hdr"]
    parts = {k: [] for k in other_keys}

    for i in tqdm(range(n_load), desc="splits"):
        split = _load_split(file, other_keys, i)
        for k in other_keys:
            if filter_df:
                split[k] = DFCleaning.filter_using_hdr(split[k], hdr_clean)
            if k == "cc1pi" and cosmic_filter:
                split[k] = _apply_obvious_cosmic_filter(split[k])   # BEFORE anything else
            if shrink:
                split[k] = shrink_df(split[k], ints=False, categories=False)
            parts[k].append(split[k])
        del split
        gc.collect()

    # ---- Concat one key at a time, then ints/categories ----
    df = {}
    for k in other_keys:
        df[k] = pd.concat(parts.pop(k))
        if shrink:
            df[k] = shrink_df(df[k], floats=False, ints=True, categories=True, verbose=True)
        gc.collect()
    if "hdr" in keys2load:
        df["hdr"] = hdr_clean

    print("loaded!")
    if reprocess_df and "cc1pi" in keys2load:
        print("Nothing to update")
    if reprocess_truth and "nudf" in keys2load:
        print("Nothing to update")
    return df
    
import gc

import pandas as pd
from tqdm import tqdm

# Assumed to exist in your project, as in your current code:
#   splh, DFCleaning, get_muon_pion_slice_mask
#   truth_cols_to_keep_slim, reco_cols_to_keep_for_data_mc_comparison,
#   min_truth_cols_to_keep, reco_cols_to_keep, min_reco_cols_to_keep,
#   reco_cols_to_keep_for_BDT_comparison

IS_PERFECT_COL = ('slc', 'cut', 'is_perfect_cc1pi', '', '', '')


def _load_split(file, keys, split_idx):
    """Load ONE split of `file` for `keys` and return a dict {key: DataFrame}.

    Mirrors splh.load_dfs, which reads the HDF key f"{key}_{i}" for each split i
    and concatenates them with the index unchanged, so concatenating the splits
    gives the same result as before.
    """
    return {k: pd.read_hdf(file, key=f"{k}_{split_idx}") for k in keys}


def slim_dfs(df, slim_key=None):
    """Slim one split (dict of DataFrames). slim_key=None -> returned untouched.

    Same rules as before, but without mutating the global column lists
    (the old code appended to them on every call).
    """
    if slim_key is None:
        return df

    has_nu = "nudf" in df
    has_reco = "cc1pi" in df
    reco_keep = list(reco_cols_to_keep)
    min_reco_keep = list(min_reco_cols_to_keep)

    if slim_key == "_for_data_mc_comp":
        if has_nu:
            df['nudf'] = df['nudf'][truth_cols_to_keep_slim]
        if has_reco:
            df['cc1pi'] = df['cc1pi'][reco_cols_to_keep_for_data_mc_comparison]

    if "_min_reco_no_syst" in slim_key and has_reco:
        df['cc1pi'][IS_PERFECT_COL] = get_muon_pion_slice_mask(df['cc1pi'])
        if IS_PERFECT_COL not in reco_keep:
            reco_keep.append(IS_PERFECT_COL)
        if IS_PERFECT_COL not in min_reco_keep:
            min_reco_keep.append(IS_PERFECT_COL)

    if "no_syst" in slim_key:
        if has_nu:
            df['nudf'] = df['nudf'][min_truth_cols_to_keep]
        if has_reco:
            df['cc1pi'] = df['cc1pi'][reco_keep]

    if "min_reco" in slim_key and has_reco:
        df['cc1pi'] = df['cc1pi'][min_reco_keep]
        df['cc1pi'] = (
            df['cc1pi']
            .groupby(['__ntuple', 'entry', 'rec.slc..index'])
            .first()
        )

    if "for_BDT_comp" in slim_key and has_reco:
        df['cc1pi'] = df['cc1pi'][reco_cols_to_keep_for_BDT_comparison]

    return df
    
# ======================================================================
# load_df_by_split
# ======================================================================
def load_df_by_split(file, keys2load, n_max_concat=100, filter_df=True,
                     reprocess_df=True, reprocess_truth=True, slim_key=None,
                     shrink=True):
    """Load a file split by split, filtering, slimming and shrinking each split as it is
    loaded so that only reduced data is ever accumulated in memory.

    slim_key=None -> no slimming (columns untouched).
    shrink=True   -> float32 per split; ints/categories after the final concat (see notes).
    """
    print(f"keys in {file}")
    splh.print_keys(file)

    n_split = splh.get_n_split(file)
    print("n_split: %d" % n_split)
    n_load = min(n_split, n_max_concat)

    print('dataframes')

    # ---- Pass 1: hdr only (small) -> global duplicate filtering ----------
    hdr_clean = None
    if filter_df:
        hdr_parts = [
            _load_split(file, ["hdr"], i)["hdr"]
            for i in tqdm(range(n_load), desc="Loading hdr", unit="split")
        ]
        hdr_all = pd.concat(hdr_parts)
        del hdr_parts
        print(hdr_all['pot'].sum())

        print(f"duplication for {file}")
        DFCleaning.find_duplicate_run_evt_combinations(hdr_all)
        DFCleaning.plot_duplicate_run_subrun_evt_distribution(hdr_all, file)

        hdr_clean = DFCleaning.filter_unique_events(hdr_all)
        DFCleaning.find_duplicate_run_evt_combinations(hdr_clean)
        del hdr_all
        gc.collect()

    # ---- Pass 2: per split -> filter with hdr -> slim -> shrink -> keep ----
    other_keys = [k for k in keys2load if k != "hdr"]
    parts = {k: [] for k in other_keys}

    for i in tqdm(range(n_load), desc="Loading + slimming splits", unit="split"):
        split = _load_split(file, other_keys, i)

        if filter_df:
            for k in other_keys:
                split[k] = DFCleaning.filter_using_hdr(split[k], hdr_clean)

        split = slim_dfs(split, slim_key)

        if shrink:
            # Per split only value-independent conversions (float64 -> float32):
            # integer widths and category sets would differ between splits and
            # pd.concat would upcast them again (categoricals with different
            # categories even fall back to object).
            for k in other_keys:
                split[k] = shrink_df(split[k], ints=False, categories=False)

        for k in other_keys:
            parts[k].append(split[k])
        del split
        gc.collect()

    # ---- Concatenate the reduced parts, one key at a time ------------------
    df = {}
    for k in other_keys:
        df[k] = pd.concat(parts.pop(k))
        if shrink:
            # Now the whole column is known: safe to pick integer widths / categories
            df[k] = shrink_df(df[k], floats=False, ints=True, categories=True, verbose=True)
        gc.collect()

    if "hdr" in keys2load:
        if hdr_clean is not None:
            df["hdr"] = hdr_clean
        else:
            df["hdr"] = pd.concat(
                [_load_split(file, ["hdr"], i)["hdr"] for i in range(n_load)]
            )
            print(df["hdr"]['pot'].sum())
    print('loaded!')

    if reprocess_df:
        if "cc1pi" in keys2load:
            print("Nothing to update")
    if reprocess_truth:
        if "nudf" in keys2load:
            print("Nothing to update")

    return df
    
import gc
import pandas as pd
from tqdm import tqdm


def perform_truth_matching(mc_evt_df, mc_nu_df):
    new_columns = []
    for c in mc_nu_df.columns:
        new_columns.append(('truth',) + c + ('',) + ('',))  # prepend 'truth'
    mc_nu_df.columns = pd.MultiIndex.from_tuples(new_columns)

    matchdf = ph.multicol_merge(
        mc_evt_df.reset_index(),
        mc_nu_df.reset_index(),
        left_on=[
            ("__ntuple", "", "", "", "", ""),
            ("entry", "", "", "", "", ""),
            ("slc", "tmatch", "idx", "", "", "")
        ],
        right_on=[
            ("__ntuple", "", "", "", "", ""),
            ("entry", "", "", "", "", ""),
            ("rec.mc.nu..index", "", "", "", "", "")
        ],
        how="left"
    )  # -- save all slices

    # Reindex so it is again "__ntuple", "entry", "slice_id"
    matchdf = matchdf.set_index(mc_evt_df.index.names, verify_integrity=True)
    
    # Remove "rec.mc.nu..index"
    matchdf = matchdf.drop(columns=[('rec.mc.nu..index', '', '', '', '', '')])

    # Helper to safely fill NaN on both Categorical and object/string series
    def _fill_cosmic(series):
        if isinstance(series.dtype, pd.CategoricalDtype):
            if 'cosmic' not in series.cat.categories:
                series = series.cat.add_categories('cosmic')
        return series.fillna('cosmic')

    target_cols = [
        ('truth', 'nu_categ', '', '', '', ''),
        ('truth', 'genie_categ', '', '', '', ''),
        ('truth', 'nu_categ_proton_reduced', '', '', '', '')
    ]

    for col in target_cols:
        if col in matchdf.columns:
            matchdf[col] = _fill_cosmic(matchdf[col])

    return matchdf

import gc
import pandas as pd
from tqdm import tqdm

def perform_truth_matching_low_memmory(mc_evt_df, mc_nu_df, ntuple_chunk_size=10):
    mc_evt_df = mc_evt_df.sort_index()
    mc_nu_df = mc_nu_df.sort_index()
    
    unique_ntuples = mc_evt_df.index.get_level_values("__ntuple").unique()
    num_ntuples = len(unique_ntuples)

    # 1. We use a generator function to yield batches
    # This prevents the creation of a massive intermediate list object
    def batch_generator():
        pbar = tqdm(range(0, num_ntuples, ntuple_chunk_size), desc="Matching Batches")
        for i in pbar:
            batch_ntuples = unique_ntuples[i : i + ntuple_chunk_size]
            
            # Use .loc with a slice if the index is sorted (much faster than .isin)
            start, end = batch_ntuples[0], batch_ntuples[-1]
            evt_chunk = mc_evt_df.loc[start:end].copy()
            nu_chunk = mc_nu_df.loc[start:end].copy()
            
            matched_batch = perform_truth_matching(evt_chunk, nu_chunk)
            
            if not matched_batch.empty:
                yield matched_batch

            # Explicitly clear chunk memory before the next yield
            del evt_chunk
            del nu_chunk
            gc.collect()

    print("Starting Truth Matching...")
    
    # 2. Use pd.concat directly on the generator
    # This is more memory efficient than building a list first
    final_df = pd.concat(batch_generator(), copy=False)
    
    # 3. Final cleanup of the generator and intermediate data
    gc.collect()
    
    print("Finalizing Result...")
    return final_df


import pandas as pd

import gc
import pandas as pd

def concat_shift_first_index(dfs, level=0, free_inputs=True):
    """
    Concatenates MultiIndex DataFrames with shifted index level, popping inputs
    sequentially to minimize peak memory overhead during concatenation.
    """
    if isinstance(level, str):
        level = dfs[0].index.names.index(level)

    parts = []
    offset = 0

    # Pop frames one-by-one from the list so original references are cleared
    while dfs:
        df = dfs.pop(0)
        if df is None or df.empty:
            continue

        lv = df.index.levels[level]
        max_val = int(lv.max())

        if offset > 0:
            # Create a shallow copy with re-leveled index (no data duplication)
            df = df.copy(deep=False)
            df.index = df.index.set_levels(lv + offset, level=level, verify_integrity=False)

        parts.append(df)
        offset += max_val + 1

    if not parts:
        return pd.DataFrame()

    # Perform final concatenation
    out = pd.concat(parts, copy=False)
    
    # Immediately clear temporary list
    parts.clear()
    gc.collect()
    
    return out

import gc
import os
import shutil
import tempfile
from concurrent.futures import (
    ProcessPoolExecutor, ThreadPoolExecutor, wait, FIRST_COMPLETED,
)

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from tqdm import tqdm

def shrink_df(df, floats=True, ints=True, categories=True, skip_cols=None,
              cat_max_ratio=0.5, verbose=False):
    """In place, column by column. The index is never modified.
       float64 -> float32 (if safe)           int64 -> int32 (if it fits)
       object of bools   -> bool / boolean    object of ints   -> int32 / Int32 (int64 if too big)
       object of floats  -> float32 (if safe) object of strings -> category (if few distinct)
    Do NOT use on the hdr DataFrame."""
    if df is None or df.empty:
        return df

    skip = set(SHRINK_SKIP_COLS) | set(skip_cols or ())
    before = df.memory_usage(deep=True).sum() if verbose else 0
    i32 = np.iinfo(np.int32)

    for c in list(df.columns):
        if _col_tokens(c) & skip:
            continue
        s = df[c]
        if not isinstance(s, pd.Series):          # duplicated label -> leave alone
            continue
        dt = s.dtype

        if floats and dt == np.float64:
            if _float_is_safe(s):
                df[c] = s.astype(np.float32)

        elif ints and pd.api.types.is_integer_dtype(dt) and dt != np.bool_:
            if np.dtype(dt).itemsize > 4 and len(s) > 0 \
               and s.min() >= i32.min and s.max() <= i32.max:
                df[c] = s.astype(np.int32)

        elif dt == object and len(s) > 0:
            kind = pd.api.types.infer_dtype(s, skipna=True)
            has_na = bool(s.isna().any())

            if kind == "boolean" and ints:
                df[c] = s.astype("boolean" if has_na else bool)

            elif kind == "integer" and ints:
                fits = s.min() >= i32.min and s.max() <= i32.max
                if has_na:
                    df[c] = s.astype("Int32" if fits else "Int64")
                else:
                    df[c] = s.astype(np.int32 if fits else np.int64)

            elif kind in ("floating", "mixed-integer-float") and floats:
                f = s.astype(np.float64)
                df[c] = f.astype(np.float32) if _float_is_safe(f) else f

            elif kind == "string" and categories \
                 and s.nunique(dropna=True) < cat_max_ratio * len(s):
                df[c] = s.astype("category")

    if verbose:
        after = df.memory_usage(deep=True).sum()
        print(f"    shrink_df: {before / 1e9:.3f} GB -> {after / 1e9:.3f} GB")
    return df
    
def _n_workers(fraction):
    try:
        n_cpus = len(os.sched_getaffinity(0))
    except AttributeError:
        n_cpus = os.cpu_count() or 1
    return max(1, int(n_cpus * fraction))


def _consolidate(paths, dest):
    # Stream chunk files into ONE parquet file, a single chunk in memory at a time.
    writer = None
    try:
        for p in tqdm(paths, desc="Consolidating", unit="file"):
            t = pq.read_table(p)
            if writer is None:
                writer = pq.ParquetWriter(dest, t.schema)
            writer.write_table(t.cast(writer.schema))
            del t
    finally:
        if writer is not None:
            writer.close()
    return dest


def _match_chunk(evt_chunk, nu_chunk, out_path):
    matched = perform_truth_matching(evt_chunk, nu_chunk)
    if matched is None or matched.empty:
        return None
    for c in matched.columns:
        if c[0] == "truth" and matched[c].dtype == np.float64:
            matched[c] = matched[c].astype(np.float32)   # same rule for every chunk -> consistent schema
    matched.to_parquet(out_path)
    return out_path
    
import os
import sys
import gc
import tempfile
import shutil
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import pyarrow.dataset as ds
import pandas as pd
from tqdm import tqdm


def perform_truth_matching_parallel(
    mc_evt_df,
    mc_nu_df,
    ntuple_chunk_size=10,
    cpu_fraction=0.8,
    max_workers=16,
    output_dir=None,
    return_df=True,
    max_tasks_per_child=3,
    consolidate_to=None,
    read_threads=8,
):
    if not return_df and consolidate_to is None and output_dir is None:
        raise ValueError("return_df=False needs output_dir or consolidate_to.")

    # Avoid making unnecessary copies if the data is already sorted
    if not mc_evt_df.index.is_monotonic_increasing:
        mc_evt_df = mc_evt_df.sort_index()
    if not mc_nu_df.index.is_monotonic_increasing:
        mc_nu_df = mc_nu_df.sort_index()

    unique_ntuples = mc_evt_df.index.get_level_values("__ntuple").unique()
    starts = list(range(0, len(unique_ntuples), ntuple_chunk_size))

    n_workers = _n_workers(cpu_fraction)
    if max_workers is not None:
        n_workers = max(1, min(n_workers, max_workers))

    print(
        f"Starting Truth Matching with {n_workers} workers, {len(starts)} chunks..."
    )

    keep_files = output_dir is not None
    out_dir = output_dir or tempfile.mkdtemp(prefix="truth_match_")
    os.makedirs(out_dir, exist_ok=True)

    def make_chunks():
        for k, i in enumerate(starts):
            batch = unique_ntuples[i : i + ntuple_chunk_size]
            start, end = batch[0], batch[-1]
            yield (
                k,
                mc_evt_df.loc[start:end],
                mc_nu_df.loc[start:end],
                os.path.join(out_dir, f"chunk_{k:06d}.parquet"),
            )

    paths = {}
    pending = {}
    chunks = make_chunks()
    exhausted = False
    pool_kwargs = {}

    # max_tasks_per_child is only supported in Python 3.11+
    if max_tasks_per_child is not None and sys.version_info >= (3, 11):
        pool_kwargs["max_tasks_per_child"] = max_tasks_per_child

    try:  # Indented at 4 spaces
        with (
            ProcessPoolExecutor(max_workers=n_workers, **pool_kwargs) as ex,
            tqdm(total=len(starts), desc="Matching Batches") as pbar,
        ):
            while pending or not exhausted:
                while not exhausted and len(pending) < n_workers:
                    try:
                        k, evt_chunk, nu_chunk, out_path = next(chunks)
                    except StopIteration:
                        exhausted = True
                        break
                    pending[
                        ex.submit(_match_chunk, evt_chunk, nu_chunk, out_path)
                    ] = k
                    del evt_chunk, nu_chunk

                if not pending:
                    break

                done, _ = wait(pending, return_when=FIRST_COMPLETED)
                for fut in done:
                    k = pending.pop(fut)
                    p = fut.result()
                    if p is not None:
                        paths[k] = p
                    pbar.update(1)

        ordered = [paths[k] for k in sorted(paths)]

        del mc_evt_df, mc_nu_df, unique_ntuples, chunks, make_chunks
        gc.collect()

        print("Finalizing Result...")
        if not ordered:
            return pd.DataFrame()

        if consolidate_to is not None:
            _consolidate(ordered, consolidate_to)
            if not return_df:
                return consolidate_to

        if not return_df:
            return ordered

        # Stream parquet dataset directly from disk
        print("Streaming Parquet dataset from disk...")
                # Read chunk files one by one with a progress bar
        tables = []
        for p in tqdm(ordered, desc="Reading chunks", unit="file"):
            tables.append(pq.read_table(p, memory_map=True))

        with tqdm(total=3, desc="Finalizing", unit="stage") as stage:
            stage.set_postfix_str("concatenating")
            table = pa.concat_tables(tables, promote_options="permissive")
            del tables
            gc.collect()
            stage.update(1)

            stage.set_postfix_str("unifying dictionaries")
            table = table.unify_dictionaries()
            stage.update(1)

            stage.set_postfix_str(
                f"to pandas ({table.num_rows:,} rows x {table.num_columns} cols)"
            )
            final_df = table.to_pandas(self_destruct=True, split_blocks=True)
            del table
            gc.collect()
            stage.update(1)

        return final_df

        print(
            f"Converting {table.num_rows:,} rows x {table.num_columns} cols to pandas..."
        )
        final_df = table.to_pandas(self_destruct=True, split_blocks=True)
        del table
        gc.collect()

        return final_df

    finally:  # Must match the 4-space indentation of 'try:'
        if not keep_files:
            shutil.rmtree(out_dir, ignore_errors=True)