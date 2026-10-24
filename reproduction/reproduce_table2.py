"""Recompute, from the printed totals of Kastrinos et al. 2021 (Gastroenterology, PMID 33839100), Table 2, the arithmetic its conclusions rest on,
with the Python standard library only.

WHAT THE TABLE IS. It compares colonoscopy schedules among carriers who are already known (a cohort of healthy 25-year-olds followed to 75, US
health-care system perspective, 3% discounting as the paper states it, 2019 US dollars, perfect adherence). It has no "not identified" arm, so nothing here is the value
of identifying a carrier. What this script reproduces is the ARITHMETIC OF THE PRINTED TOTALS; it does not rebuild the paper's Markov model,
whose parameters are not in the table.

1. Every printed strategy is classified: on the efficiency frontier, or dominated (another strategy costs no more and gains no fewer QALYs).
2. Each printed ICER is set beside the one recomputed from the printed totals, allowing for the rounding of QALYs to three decimals.
3. The optimal strategy at 100,000 USD per QALY (the paper's own definition: the highest QALYs with an ICER no greater than the threshold), and its
   net monetary benefit, NMB = threshold x incremental QALYs - incremental cost, against two comparators that answer different questions: the next
   efficient strategy on the frontier, and the least intensive strategy listed (current management for a family history, which the paper rates dominated).
4. A check of the paper's statement that costs and QALYs are discounted at 3% a year: the most 3% discounting can give a cohort followed for 50 annual
   cycles (ages 25 to 75), at full health with no deaths, is compared with the printed QALYs.

Usage: python3 reproduction/reproduce_table2.py   (the output is compared with EXPECTED_OUTPUT.txt by tests/public)."""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WTP = 100_000.0
MISSING = ("kastrinos2021_table2.tsv is not shipped with this repository: the publisher's terms for the NIH author manuscript do not grant redistribution.\n"
           "See reproduction/README.md for how to supply it (five columns, 20 printed strategies). Nothing was computed.\n")


def money(s):
    return float(str(s).replace("$", "").replace(",", ""))


def table():
    if not os.path.exists(os.path.join(HERE, "kastrinos2021_table2.tsv")):
        sys.stderr.write(MISSING)
        raise SystemExit(2)
    lines = [l for l in open(os.path.join(HERE, "kastrinos2021_table2.tsv"), encoding="utf-8").read().splitlines() if l and not l.startswith("#")]
    out = {}
    for r in csv.DictReader(lines, delimiter="\t"):
        out.setdefault(r["gene"], []).append({"name": r["strategy_as_printed"], "qaly": float(r["QALYs"]), "cost": money(r["Cost"]), "icer": r["ICER"]})
    return out


def classify(rows):
    n = len(rows)
    strong = [any(rows[j]["cost"] <= rows[i]["cost"] and rows[j]["qaly"] >= rows[i]["qaly"]
                  and (rows[j]["cost"] < rows[i]["cost"] or rows[j]["qaly"] > rows[i]["qaly"]) for j in range(n) if j != i) for i in range(n)]
    alive = sorted((i for i in range(n) if not strong[i]), key=lambda i: rows[i]["cost"])
    ext = set()
    while True:
        for k in range(1, len(alive) - 1):
            a, b, c = rows[alive[k - 1]], rows[alive[k]], rows[alive[k + 1]]
            if (b["cost"] - a["cost"]) / (b["qaly"] - a["qaly"]) > (c["cost"] - b["cost"]) / (c["qaly"] - b["qaly"]):
                ext.add(alive[k])
                alive.pop(k)
                break
        else:
            break
    return ["dominated" if strong[i] else "extended-dominated" if i in ext else "frontier" for i in range(n)], alive


def sgn(x, nd=2):
    return ("-" if x < 0 else "+") + f"{abs(x):,.{nd}f}"


def main():
    res = {g: (rows, *classify(rows)) for g, rows in table().items()}
    print("1. Classification of every printed strategy")
    for g, (rows, status, front) in res.items():
        for i, r in enumerate(rows):
            print(f"  {g} {r['name']}: {status[i]} (printed {r['icer']})")
    print("2. Printed ICERs against the previous strategy on the frontier, recomputed from the printed totals")
    steps = {}
    for g, (rows, status, front) in res.items():
        for k in range(1, len(front)):
            a, b = rows[front[k - 1]], rows[front[k]]
            dq, dc = b["qaly"] - a["qaly"], b["cost"] - a["cost"]
            lo, hi = dc / (dq + 0.001), dc / (dq - 0.001)
            printed = money(b["icer"])
            ok = lo - 1 <= printed <= hi + 1
            steps[(g, k)] = dc / dq
            line = f"  {g} {b['name']} vs {a['name']}: printed {printed:,.2f}, recomputed {dc / dq:,.2f}, "
            line += "reproduces" if ok else (f"does not reproduce (QALYs rounded to three decimals allow {lo:,.2f} to {hi:,.2f}; "
                                              f"the printed ratio implies an increment of {dc / printed:.4f} QALYs, the printed totals give {dq:.4f})")
            print(line)
    print("3. The optimal strategy at 100,000 USD per QALY, and its net monetary benefit against two comparators")
    for g, (rows, status, front) in res.items():
        elig = [front[0]] + [front[k] for k in range(1, len(front)) if steps[(g, k)] <= WTP]
        opt = max(elig, key=lambda i: rows[i]["qaly"])
        o = rows[opt]
        print(f"  {g} optimal: {o['name']}")
        k = front.index(opt)
        if k > 0:
            p = rows[front[k - 1]]
            dq, dc = o["qaly"] - p["qaly"], o["cost"] - p["cost"]
            print(f"    against the next efficient strategy ({p['name']}): QALYs {sgn(dq, 3)}, cost {sgn(dc)}, NMB {WTP * dq - dc:,.2f}")
        else:
            print("    against the next efficient strategy: none (the optimal strategy is the reference of the frontier)")
        ref = rows[0]
        dq, dc = o["qaly"] - ref["qaly"], o["cost"] - ref["cost"]
        print(f"    against the least intensive strategy listed ({ref['name']}): QALYs {sgn(dq, 3)}, cost {sgn(dc)}, NMB {WTP * dq - dc:,.2f}")
    print("4. A check of the stated discounting (the source states 3% a year, ages 25 to 75)")
    cap = sum(1.03 ** -t for t in range(50))
    cap_end = sum(1.03 ** -(t + 1) for t in range(50))
    print(f"  The most 3% discounting can give over 50 annual cycles, at full health with no deaths, is {cap:.2f} QALYs ({cap_end:.2f} if the first cycle is also discounted)")
    allover = True
    for g, (rows, status, front) in res.items():
        qs = [r["qaly"] for r in rows]
        allover = allover and min(qs) > cap
        print(f"  {g}: printed QALYs {min(qs):.3f} to {max(qs):.3f}, {'all above' if min(qs) > cap else 'not all above'} {cap:.2f}")
    if allover:
        print("  The printed QALYs therefore cannot be 3% discounted over that horizon: the discounting convention behind the printed totals is unresolved, and every figure above is arithmetic on the totals as printed.")


if __name__ == "__main__":
    main()
