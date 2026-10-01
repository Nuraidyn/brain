"""Парсер readable-логов PsychoPy -> файлы таймингов для SPM25 (fMRI Model Specification).

Использование:
    python parser.py <input_readable.txt> [--design lang|switch] [--pmod-rt]
                     [--out DIR] [--decimals N]

Для каждой сессии (RUN1 = Block1_ColorWord, RUN2 = Block2_ConceptWord) создаётся
MATLAB-скрипт <participant>_event_SPM_<design>_run<N>.txt с условиями одного дизайна:
    onsets = {[...], ...};
    durations = {[...], ...};
    names = {'...', ...};
    save <participant>_event_spm_<design>_run<N>.mat names onsets durations;
"""

import argparse
import ast
import re
import sys
from pathlib import Path

SESSIONS = {"Block1_ColorWord": 1, "Block2_ConceptWord": 2}

ONSET_KEYS = ("MixS.started", "mixC.started", "thisRow.t")
STOP_KEYS = ("MixS.stopped", "mixC.stopped")

LANG_ABBR = {"English": "Eng", "Kazakh": "Kaz", "Russian": "Rus"}
LANG_CONDITIONS = [f"{l}_{c}" for l in ("Eng", "Kaz", "Rus") for c in ("Congruent", "Incongruent")]
SWITCH_CONDITIONS = ["Switch", "Repeat", "Filler"]

RECORD_RE = re.compile(r"^RECORD (\d+): (.+)$")


# ---------------------------------------------------------------- чтение файла

def parse_file(path):
    """Возвращает (metadata, records), где records — список (номер, заголовок, dict)."""
    lines = Path(path).read_text(encoding="utf-8-sig").splitlines()

    metadata = {}
    records = []
    current_header = None
    buffer = []

    def flush():
        if current_header is None or not buffer:
            return
        body = ast.literal_eval("\n".join(buffer))
        if current_header == "METADATA":
            metadata.update(body)
        else:
            records.append((*current_header, body))

    for line in lines:
        m = RECORD_RE.match(line)
        if m or line.strip() == "EXPERIMENT METADATA":
            flush()
            buffer = []
            current_header = (int(m.group(1)), m.group(2)) if m else "METADATA"
        elif current_header is not None and line.strip():
            buffer.append(line)
    flush()
    return metadata, records


def first_value(record, keys):
    for k in keys:
        v = record.get(k)
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            return float(v)
    return None


# ---------------------------------------------------------------- логика

def collect_sessions(records, warnings):
    """Делит экспериментальные триалы по сессиям и пересчитывает onset от InitialFix_RUNx."""
    t0 = {}
    trials = {run: [] for run in SESSIONS.values()}

    for num, _, rec in records:
        for run in SESSIONS.values():
            key = f"InitialFix_RUN{run}.started"
            if key in rec and run not in t0:
                t0[run] = float(rec[key])
        run = SESSIONS.get(rec.get("Block"))
        if run is None:
            continue  # PB_S, PB_C, Other routine
        onset = first_value(rec, ONSET_KEYS)
        if onset is None:
            warnings.append(f"RECORD {num}: нет времени onset, триал пропущен")
            continue
        stop = first_value(rec, STOP_KEYS)
        started = first_value(rec, ONSET_KEYS[:2])
        if stop is None or started is None:
            warnings.append(f"RECORD {num}: нет stopped/started, duration = 0")
            duration = 0.0
        else:
            duration = stop - started
        trials[run].append({
            "record": num,
            "onset": onset,
            "duration": duration,
            "language": rec.get("Language"),
            "congruent": rec.get("congruent"),
            "transition": rec.get("transition"),
            "correct": rec.get("luminaCorrect"),
            "rt": rec.get("luminaRT"),
        })

    for run, items in trials.items():
        if not items:
            continue
        if run not in t0:
            t0[run] = items[0]["onset"]
            warnings.append(f"Сессия {run}: InitialFix_RUN{run} не найден, t0 = onset первого триала")
        for t in items:
            t["onset"] -= t0[run]
    return trials


def condition_of(trial, design):
    if design == "lang":
        lang = LANG_ABBR.get(trial["language"])
        if lang is None or trial["congruent"] not in (0, 1):
            return None
        return f"{lang}_{'Congruent' if trial['congruent'] == 1 else 'Incongruent'}"
    transition = trial["transition"]
    return transition.capitalize() if isinstance(transition, str) else None


def build_design(trials, design, run, warnings):
    order = LANG_CONDITIONS if design == "lang" else SWITCH_CONDITIONS
    groups = {name: [] for name in order}
    for t in trials:
        if t["correct"] != 1:
            continue  # ошибки и пропуски исключаются
        name = condition_of(t, design)
        if name not in groups:
            warnings.append(f"RECORD {t['record']}: не удалось определить условие ({design})")
            continue
        groups[name].append(t)
    for name in order:
        if not groups[name]:
            warnings.append(f"Сессия {run}, дизайн {design}: условие {name} пустое, пропущено")
    return {name: items for name, items in groups.items() if items}


# ---------------------------------------------------------------- вывод

def fmt_vector(values, decimals):
    return "[" + ", ".join(f"{v:.{decimals}f}" for v in values) + "]"


def write_matlab(path, mat_name, groups, decimals, pmod_rt):
    names = list(groups)
    onsets = [fmt_vector([t["onset"] for t in groups[n]], decimals) for n in names]
    durations = [fmt_vector([t["duration"] for t in groups[n]], decimals) for n in names]

    lines = [
        "onsets = {" + ", ".join(onsets) + "};",
        "durations = {" + ", ".join(durations) + "};",
        "names = {" + ", ".join(f"'{n}'" for n in names) + "};",
    ]
    save_vars = "names onsets durations"
    if pmod_rt:
        lines.append("pmod = struct('name', {}, 'param', {}, 'poly', {});")
        for i, n in enumerate(names, 1):
            rts = [t["rt"] if isinstance(t["rt"], (int, float)) else 0.0 for t in groups[n]]
            lines.append(f"pmod({i}).name{{1}} = 'RT'; pmod({i}).param{{1}} = "
                         f"{fmt_vector(rts, 4)}; pmod({i}).poly{{1}} = 1;")
        save_vars += " pmod"
    lines.append(f"save {mat_name} {save_vars};")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def try_write_mat(path, groups, pmod_rt):
    """Дополнительно пишет .mat напрямую, если установлен scipy."""
    try:
        import numpy as np
        from scipy.io import savemat
    except ImportError:
        return False
    names = list(groups)
    obj = lambda xs: np.array(xs, dtype=object)
    data = {
        "names": obj(names),
        "onsets": obj([np.array([t["onset"] for t in groups[n]]) for n in names]),
        "durations": obj([np.array([t["duration"] for t in groups[n]]) for n in names]),
    }
    if pmod_rt:
        data["pmod"] = np.array(
            [(obj(["RT"]),
              obj([np.array([t["rt"] if isinstance(t["rt"], (int, float)) else 0.0
                             for t in groups[n]])]),
              obj([1]))
             for n in names],
            dtype=[("name", object), ("param", object), ("poly", object)])
    savemat(path, data)
    return True


# ---------------------------------------------------------------- main

def main(argv=None):
    ap = argparse.ArgumentParser(description="PsychoPy readable log -> SPM25 timings")
    ap.add_argument("input", help="файл *_readable.txt")
    # дизайны не смешиваются: оба покрывают одни и те же триалы -> вырожденная матрица
    ap.add_argument("--design", choices=("lang", "switch"), default="lang",
                    help="lang = язык x конгруэнтность (по умолчанию), switch = Switch/Repeat/Filler")
    ap.add_argument("--pmod-rt", action="store_true",
                    help="добавить luminaRT как параметрический модулятор")
    ap.add_argument("--out", help="папка для результатов (по умолчанию рядом с входным файлом)")
    ap.add_argument("--decimals", type=int, default=1, help="знаков после запятой (по умолчанию 1)")
    args = ap.parse_args(argv)

    in_path = Path(args.input)
    out_dir = Path(args.out) if args.out else in_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    metadata, records = parse_file(in_path)
    participant = str(metadata.get("participant") or in_path.name.split("_")[0])
    warnings = []
    sessions = collect_sessions(records, warnings)

    print(f"Испытуемый: {participant} | expName: {metadata.get('expName')} | "
          f"frameRate: {metadata.get('frameRate')}")

    summary = []
    for run, trials in sessions.items():
        if not trials:
            warnings.append(f"Сессия {run}: триалы не найдены")
            continue
        errors = sum(1 for t in trials if t["correct"] == 0)
        misses = sum(1 for t in trials if t["correct"] is None)
        included = sum(1 for t in trials if t["correct"] == 1)
        summary.append(f"Сессия {run} — {len(trials)} триалов "
                       f"(включено {included}, ошибок {errors}, пропусков {misses})")
        groups = build_design(trials, args.design, run, warnings)
        base = f"{participant}_event_SPM_{args.design}_run{run}"
        mat_name = f"{participant}_event_spm_{args.design}_run{run}.mat"
        write_matlab(out_dir / f"{base}.txt", mat_name, groups, args.decimals, args.pmod_rt)
        mat_written = try_write_mat(out_dir / mat_name, groups, args.pmod_rt)
        counts = ", ".join(f"{n}: {len(g)}" for n, g in groups.items())
        print(f"  RUN{run} [{args.design}] -> {base}.txt"
              f"{' + ' + mat_name if mat_written else ''} ({counts})")

    print(f"Успешно обработан испытуемый {participant}: " + ", ".join(summary) + ".")
    if warnings:
        print(f"Предупреждения ({len(warnings)}):")
        for w in warnings:
            print("  - " + w)
    else:
        print("Ошибок не обнаружено.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
