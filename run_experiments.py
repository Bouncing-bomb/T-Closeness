#!/usr/bin/env python3
"""Reproducible, read-only Adult experiments and dependency-free SVG figures."""

import argparse
import contextlib
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import random
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
import unittest
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
import mondrian_tcloseness as algorithm
from utils.read_adult_data import ATT_NAMES, QI_INDEX, read_data

ROOT = Path(__file__).resolve().parent
SEED = 20261102
SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)
METRICS = [
    ("gcp_percent", "GCP (%)"),
    ("max_emd", "Maximum EMD"),
    ("runtime_seconds", "Partition runtime (s)"),
    ("released_ec_count", "Published equivalence class count"),
]


class AlgorithmChecks(unittest.TestCase):
    """Small cases exercise the supplied TODO contracts without writing data."""

    def setUp(self):
        self.saved_t = algorithm.T_THRESHOLD

    def tearDown(self):
        algorithm.T_THRESHOLD = self.saved_t

    def partition(self, data, k=2, qi=-1):
        algorithm.init(data, k, qi)
        return algorithm.Partition(
            data, [0] * algorithm.QI_LEN,
            [len(order) - 1 for order in algorithm.QI_ORDER],
        )

    def split(self, data, k, t, relaxed):
        partition = self.partition(data, k)
        algorithm.T_THRESHOLD = t
        method = algorithm.anonymize_relaxed if relaxed else algorithm.anonymize_strict
        method(partition)
        members = [record for group in algorithm.RESULT for record in group.member]
        self.assertEqual(Counter(map(tuple, members)), Counter(map(tuple, data)))
        self.assertEqual(len({id(record) for record in members}), len(data))
        for group in algorithm.RESULT:
            self.assertGreaterEqual(len(group), k)
            self.assertLessEqual(algorithm.compute_emd(group), t + 1e-9)
            for record in group.member:
                for dim in range(algorithm.QI_LEN):
                    index = algorithm.QI_DICT[dim][record[dim]]
                    self.assertLessEqual(group.low[dim], index)
                    self.assertLessEqual(index, group.high[dim])
        return algorithm.RESULT

    def test_global_distribution_and_reset(self):
        self.partition([[i, sa] for i, sa in enumerate("AABCCC")])
        self.assertEqual(algorithm.GLOBAL_SA_DIST, {"A": 2 / 6, "B": 1 / 6, "C": 3 / 6})
        self.partition([[0, "Z"], [1, "Z"]])
        self.assertEqual(algorithm.GLOBAL_SA_DIST, {"Z": 1.0})

    def test_emd_full_missing_and_empty_support(self):
        whole = self.partition([[i, sa] for i, sa in enumerate("AABCCC")])
        self.assertAlmostEqual(algorithm.compute_emd(whole), 0.0)
        single = algorithm.Partition([[2, "B"]], [0], [5])
        self.assertAlmostEqual(algorithm.compute_emd(single), 5 / 6)
        absent = algorithm.Partition([[2, "D"]], [0], [5])
        self.assertAlmostEqual(algorithm.compute_emd(absent), 1.0)
        empty = algorithm.Partition([], [0], [5])
        self.assertEqual(algorithm.compute_emd(empty), 0.0)
        self.assertTrue(algorithm.satisfy_t_closeness(empty, 0))

    def test_emd_threshold_equality(self):
        self.partition([[0, "A"], [1, "B"]])
        group = algorithm.Partition([[0, "A"]], [0], [0])
        self.assertTrue(algorithm.satisfy_t_closeness(group, 0.5))
        self.assertFalse(algorithm.satisfy_t_closeness(group, 0.49))

    def test_reject_both_modes(self):
        for relaxed in (False, True):
            data = [[i, "A" if i < 4 else "B"] for i in range(8)]
            self.assertEqual(len(self.split(data, 2, 0.0, relaxed)), 1)

    def test_retry_another_dimension_both_modes(self):
        for relaxed in (False, True):
            data = [[i, (i % 4) * 2 + i // 4, "A" if i < 4 else "B"]
                    for i in range(8)]
            groups = self.split(data, 2, 0.0, relaxed)
            self.assertEqual(len(groups), 4)

    def test_strict_ties_and_relaxed_balance(self):
        values = [1, 2, 3, 3, 3, 4, 5]
        strict = self.split([[v, "A"] for v in values], 3, 1.0, False)
        self.assertEqual([len(group) for group in strict], [7])
        relaxed = self.split([[v, "A"] for v in values], 3, 1.0, True)
        self.assertEqual(sorted(len(group) for group in relaxed), [3, 4])

    def test_odd_size_and_unsplittable_partition(self):
        for relaxed in (False, True):
            groups = self.split([[i, "A"] for i in range(5)], 2, 1.0, relaxed)
            self.assertEqual(sorted(len(group) for group in groups), [2, 3])
            groups = self.split([[i, "A"] for i in range(3)], 2, 1.0, relaxed)
            self.assertEqual(len(groups), 1)

    def test_constant_qi(self):
        for relaxed in (False, True):
            groups = self.split([[1, sa] for sa in "ABAB"], 2, 1.0, relaxed)
            self.assertEqual(len(groups), 1)

    def test_randomized_invariants(self):
        rng = random.Random(SEED)
        for relaxed in (False, True):
            for size in (12, 21, 40):
                for t in (0.0, 0.15, 0.5, 1.0):
                    data = [[rng.randrange(6), rng.randrange(9), rng.choice("ABC")]
                            for _ in range(size)]
                    self.split(data, 3, t, relaxed)


def run_checks():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(AlgorithmChecks)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise RuntimeError("Algorithm checks failed; experiments were not run.")
    return result.testsRun


def distance(counts, total, reference):
    return 0.5 * sum(abs(counts.get(sa, 0) / total - reference.get(sa, 0.0))
                     for sa in set(counts) | set(reference))


def measure(source, k, qi, t, relaxed):
    # mondrian mutates QI cells: fresh row copies isolate every measurement.
    data = [record[:] for record in source]
    algorithm.T_THRESHOLD = t
    original_counts = Counter(record[-1] for record in data)
    reference = {sa: count / len(data) for sa, count in original_counts.items()}
    with contextlib.redirect_stdout(io.StringIO()):
        output, (gcp, runtime) = algorithm.mondrian(data, k, relaxed, qi)
    groups = algorithm.RESULT
    assert len(output) == len(source)
    assert Counter(record[-1] for record in output) == original_counts
    assert len({id(record) for group in groups for record in group.member}) == len(source)
    assert sum(map(len, groups)) == len(source)
    emds = [distance(Counter(record[-1] for record in group.member), len(group), reference)
            for group in groups]
    violations = sum(emd > t + 1e-9 for emd in emds)
    assert violations == 0
    assert all(len(group) >= k for group in groups)
    assert all(math.isclose(algorithm.compute_emd(group), emd, abs_tol=1e-12)
               for group, emd in zip(groups, emds))
    expected_gcp = 100 * sum(
        len(group) * sum(algorithm.get_normalized_width(group, dim) for dim in range(qi))
        for group in groups
    ) / (len(source) * qi)
    assert math.isclose(gcp, expected_gcp, abs_tol=1e-10)
    # Relaxed partitions can have identical published QI intervals; count their union too.
    released = defaultdict(Counter)
    for record in output:
        released[tuple(record[:qi])][record[-1]] += 1
    released_emds = [distance(counts, sum(counts.values()), reference)
                     for counts in released.values()]
    assert all(sum(counts.values()) >= k for counts in released.values())
    assert all(emd <= t + 1e-9 for emd in released_emds)
    return {
        "gcp_percent": gcp,
        "max_emd": max(emds, default=0.0),
        "runtime_seconds": runtime,
        "partition_count": len(groups),
        "released_ec_count": len(released),
        "released_max_emd": max(released_emds, default=0.0),
        "min_partition_size": min(map(len, groups)),
        "violations": violations,
    }


def svg_node(parent, tag, attributes=None, text=None):
    node = ET.SubElement(parent, f"{{{SVG_NS}}}{tag}", attributes or {})
    if text is not None:
        node.text = str(text)
    return node


def svg_canvas(width, height, title):
    root = ET.Element(f"{{{SVG_NS}}}svg", {
        "width": str(width), "height": str(height),
        "viewBox": f"0 0 {width} {height}", "role": "img",
    })
    svg_node(root, "title", text=title)
    svg_node(root, "rect", {"width": "100%", "height": "100%", "fill": "#f8fafc"})
    svg_node(root, "style", text="text{font-family:Arial,sans-serif;fill:#334155;font-size:12px}")
    return root


def save_svg(root, path):
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def plot_sweep(summary, sweep, output):
    xlabel = {"k": "k", "qi": "Number of QIs (prefix)", "t": "t threshold",
              "datasize": "Number of records"}[sweep]
    rows = [row for row in summary if row["sweep"] == sweep]
    xvalues = sorted({row["value"] for row in rows})
    root = svg_canvas(1100, 720, f"Mondrian t-closeness: {sweep} sweep")
    svg_node(root, "text", {"x": "38", "y": "32", "style": "font-size:21px;font-weight:bold"},
             f"Effect of {xlabel}")
    svg_node(root, "text", {"x": "38", "y": "54"},
             "Fixed seed; equally spaced tested configurations; runtime: median with min-max bars")
    for mode, color, x in (("strict", "#2563eb", 800), ("relaxed", "#e76f51", 925)):
        svg_node(root, "line", {"x1": str(x), "x2": str(x + 25), "y1": "33", "y2": "33",
                                "stroke": color, "stroke-width": "3"})
        svg_node(root, "text", {"x": str(x + 31), "y": "37"}, mode)
    for panel, (key, ylabel) in enumerate(METRICS):
        left = 80 + (panel % 2) * 540
        top = 108 + (panel // 2) * 310
        width, height = 430, 205
        maximum = max(row[key] for row in rows)
        if key == "runtime_seconds":
            maximum = max(row["runtime_max"] for row in rows)
        if sweep == "t" and key == "max_emd":
            maximum = max(maximum, max(xvalues))
        ymax = maximum * 1.12 if maximum else 1.0
        # Evenly spaced configurations make all tested thresholds readable.
        def xp(index):
            return left + index * width / max(1, len(xvalues) - 1)

        def yp(number):
            return top + height - number / ymax * height

        svg_node(root, "text", {"x": str(left), "y": str(top - 15),
                                "style": "font-size:16px;font-weight:bold"}, ylabel)
        for tick in range(5):
            val = ymax * tick / 4
            y = yp(val)
            svg_node(root, "line", {"x1": str(left), "x2": str(left + width), "y1": str(y),
                                    "y2": str(y), "stroke": "#e2e8f0"})
            label = f"{val:.3f}" if key in ("max_emd", "runtime_seconds") else f"{val:.1f}"
            svg_node(root, "text", {"x": str(left - 10), "y": str(y + 4),
                                    "text-anchor": "end"}, label)
        for index, val in enumerate(xvalues):
            svg_node(root, "text", {"x": str(xp(index)), "y": str(top + height + 22),
                                    "text-anchor": "middle"}, f"{val:g}")
        svg_node(root, "text", {"x": str(left + width / 2), "y": str(top + height + 45),
                                "text-anchor": "middle"}, xlabel)
        for mode, color in (("strict", "#2563eb"), ("relaxed", "#e76f51")):
            series = sorted((row for row in rows if row["mode"] == mode), key=lambda r: r["value"])
            points = " ".join(f"{xp(i)},{yp(row[key])}" for i, row in enumerate(series))
            svg_node(root, "polyline", {"points": points, "fill": "none", "stroke": color,
                                        "stroke-width": "2.5"})
            for index, row in enumerate(series):
                x, y = xp(index), yp(row[key])
                if key == "runtime_seconds":
                    svg_node(root, "line", {"x1": str(x), "x2": str(x),
                                            "y1": str(yp(row["runtime_min"])),
                                            "y2": str(yp(row["runtime_max"])), "stroke": color})
                svg_node(root, "circle", {"cx": str(x), "cy": str(y), "r": "3.5", "fill": color})
        if sweep == "t" and key == "max_emd":
            points = " ".join(f"{xp(i)},{yp(val)}" for i, val in enumerate(xvalues))
            svg_node(root, "polyline", {"points": points, "fill": "none", "stroke": "#64748b",
                                        "stroke-dasharray": "5,5"})
            svg_node(root, "text", {"x": str(left + 235), "y": str(top - 15)},
                     "dashed: target t")
    save_svg(root, output / f"{sweep}_effects.svg")


def plot_flowchart(output):
    root = svg_canvas(1040, 1150, "Mondrian with k-anonymity and t-closeness flowchart")
    defs = svg_node(root, "defs")
    marker = svg_node(defs, "marker", {"id": "arrow", "viewBox": "0 0 10 10", "refX": "9",
                                      "refY": "5", "markerWidth": "7", "markerHeight": "7",
                                      "orient": "auto-start-reverse"})
    svg_node(marker, "path", {"d": "M 0 0 L 10 5 L 0 10 z", "fill": "#475569"})

    def box(x, y, width, lines, decision=False):
        height = 70
        if decision:
            points = f"{x + width / 2},{y - 15} {x + width},{y + 35} {x + width / 2},{y + 85} {x},{y + 35}"
            svg_node(root, "polygon", {"points": points, "fill": "#dbeafe", "stroke": "#2563eb"})
        else:
            svg_node(root, "rect", {"x": str(x), "y": str(y), "width": str(width),
                                    "height": str(height), "rx": "10", "fill": "#eef2ff",
                                    "stroke": "#6366f1"})
        for index, line in enumerate(lines):
            svg_node(root, "text", {"x": str(x + width / 2),
                                    "y": str(y + 29 + index * 18), "text-anchor": "middle"}, line)

    def arrow(points, label=None, lx=0, ly=0):
        svg_node(root, "polyline", {"points": points, "fill": "none", "stroke": "#475569",
                                    "stroke-width": "1.8", "marker-end": "url(#arrow)"})
        if label:
            svg_node(root, "text", {"x": str(lx), "y": str(ly)}, label)

    svg_node(root, "text", {"x": "520", "y": "30", "text-anchor": "middle",
                            "style": "font-size:22px;font-weight:bold"}, "Mondrian + t-closeness")
    box(330, 50, 380, ["Input: numeric/encoded QIs, last-column SA", "Valid input: N >= k >= 1; 0 <= t <= 1"])
    box(330, 150, 380, ["Initialize QI orders, ranges and global SA Q", "Root partition = all records"])
    box(330, 260, 380, ["Any allowed QI dimension left?"], True)
    box(330, 375, 380, ["Choose allowed dimension with largest NCP", "Find median; tighten this dimension's bounds"])
    box(330, 490, 380, ["Usable median / distinct QI values?"], True)
    box(330, 610, 380, ["Strict: left <= median, right > median", "Relaxed: distribute ties to balance counts"])
    box(330, 725, 380, ["Both child sizes >= k?"], True)
    box(330, 845, 380, ["Both EMD(child, Q) <= t + 1e-9?"], True)
    box(330, 965, 380, ["Accept split; recursively anonymize both children", "New children start with all dimensions allowed"])
    box(20, 575, 235, ["Reject current cut", "Set allow[dimension] = 0"])
    box(780, 260, 240, ["Save parent as terminal partition", "Generalize QIs; output metrics"])
    arrow("520,120 520,150")
    arrow("520,220 520,245")
    arrow("520,345 520,375", "Yes", 532, 364)
    arrow("710,295 780,295", "No", 737, 285)
    arrow("520,445 520,475")
    arrow("520,575 520,610", "Yes", 532, 595)
    arrow("520,680 520,710")
    arrow("520,810 520,830", "Yes", 532, 825)
    arrow("520,930 520,965", "Yes", 532, 950)
    arrow("330,525 137,525 137,575", "No", 235, 513)
    arrow("330,760 137,760 137,645", "No", 235, 748)
    arrow("330,880 80,880 80,645", "No", 235, 868)
    arrow("137,575 137,295 330,295", "Retry same parent", 144, 370)
    arrow("710,1000 750,1000 750,235 520,235 520,245")
    svg_node(root, "text", {"x": "775", "y": "995"}, "For each child")
    svg_node(root, "text", {"x": "40", "y": "1090"},
             "Retry implementation: strict = continue; relaxed = recursive retry and return.")
    svg_node(root, "text", {"x": "40", "y": "1115"},
             "Categorical SA with unit ground distance: EMD = 0.5 * sum_v |P(v) - Q(v)|.")
    save_svg(root, output / "algorithm_flowchart.svg")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path, default=ROOT / "experiment_results")
    parser.add_argument("--checks-only", action="store_true")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    checks = run_checks()
    if args.checks_only:
        return
    output = args.output.resolve()
    for protected in (ROOT / "data", ROOT / "utils"):
        if output == protected or protected in output.parents:
            parser.error("Experiment output must not be inside data/ or utils/")
    output.mkdir(parents=True, exist_ok=True)
    original_cwd = Path.cwd()
    original_t = algorithm.T_THRESHOLD
    try:
        # The unchanged loader uses a relative input path; it only reads adult.data.
        os.chdir(ROOT)
        data, _ = read_data()
        random.Random(SEED).shuffle(data)
        size = len(data)
        qi = len(QI_INDEX)
        baseline = {"k": 10, "qi": qi, "t": 0.1, "datasize": size}
        sweeps = {
            "k": [5, 10, 20, 50, 100],
            "qi": list(range(1, qi + 1)),
            "t": [0.0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0],
            "datasize": sorted({n for n in [1000, 5000, 10000, 20000, size] if n <= size}),
        }
        raw = []
        cache = {}
        for sweep, values in sweeps.items():
            for value in values:
                config = dict(baseline, **{sweep: value})
                source = data[:config["datasize"]]
                for relaxed in (False, True):
                    mode = "relaxed" if relaxed else "strict"
                    key = (mode, config["k"], config["qi"], config["t"], config["datasize"])
                    if key not in cache:
                        # One unrecorded warmup per unique configuration.
                        measure(source, config["k"], config["qi"], config["t"], relaxed)
                        cache[key] = [measure(source, config["k"], config["qi"], config["t"], relaxed)
                                      for _ in range(args.repeats)]
                    for repeat, metrics in enumerate(cache[key], 1):
                        raw.append({"sweep": sweep, "value": value, "mode": mode,
                                    **config, "repeat": repeat, **metrics})
                    metrics = cache[key][0]
                    print(f"{sweep}={value:g} {mode}: GCP={metrics['gcp_percent']:.4f}% "
                          f"max_EMD={metrics['max_emd']:.6f} partitions={metrics['partition_count']}",
                          flush=True)
        with (output / "measurements.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(raw[0]))
            writer.writeheader()
            writer.writerows(raw)
        grouped = defaultdict(list)
        for row in raw:
            grouped[(row["sweep"], row["value"], row["mode"])].append(row)
        summary = []
        for (sweep, val, mode), measurements in grouped.items():
            row = {"sweep": sweep, "value": val, "mode": mode,
                   **{name: measurements[0][name] for name in baseline}}
            for name in cache[next(iter(cache))][0]:
                if name == "runtime_seconds":
                    times = [m[name] for m in measurements]
                    row[name] = statistics.median(times)
                    row["runtime_min"] = min(times)
                    row["runtime_max"] = max(times)
                    row["runtime_stdev"] = statistics.stdev(times) if len(times) > 1 else 0.0
                else:
                    values = [m[name] for m in measurements]
                    assert all(math.isclose(v, values[0], abs_tol=1e-12) for v in values)
                    row[name] = values[0]
            summary.append(row)
        metadata = {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "python": sys.version, "platform": platform.platform(),
            "machine": platform.machine(), "processor": platform.processor(),
            "seed": SEED, "repeats": args.repeats, "warmups_per_configuration": 1,
            "checks_passed": checks, "unique_configurations": len(cache),
            "measured_runs": len(cache) * args.repeats,
            "baseline": baseline, "sweeps": sweeps,
            "qi_names": [ATT_NAMES[i] for i in QI_INDEX], "sensitive_attribute": "income",
            "records": size, "global_sa_counts": dict(Counter(r[-1] for r in data)),
            "data_sha256": hashlib.sha256((ROOT / "data/adult.data").read_bytes()).hexdigest(),
            "algorithm_sha256": hashlib.sha256((ROOT / "mondrian_tcloseness.py").read_bytes()).hexdigest(),
            "runtime_definition": "mondrian rtime: partition recursion only; excludes init, generalization and metrics",
            "sampling": "one seeded shuffle, nested prefixes; fixed rows and order across repeats",
            "baseline_rows_reused_across_sweeps": True,
        }
        (output / "summary.json").write_text(json.dumps(
            {"metadata": metadata, "results": summary}, ensure_ascii=False, indent=2,
        ) + "\n", encoding="utf-8")
        for sweep in sweeps:
            plot_sweep(summary, sweep, output)
        plot_flowchart(output)
        print(f"Saved {len(cache)} configurations to {output}")
    finally:
        algorithm.T_THRESHOLD = original_t
        os.chdir(original_cwd)


if __name__ == "__main__":
    main()
