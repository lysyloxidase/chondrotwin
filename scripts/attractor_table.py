#!/usr/bin/env python
"""Print the attractors of the Boolean chondrocyte network for every scenario (Markdown table).

Usage:  python scripts/attractor_table.py > docs/attractors.md
"""
from chondrotwin.boolean import attractor_table

if __name__ == "__main__":
    rows = attractor_table()
    print("| Scenario | Attractor | Fate label | Basin fraction | Active nodes |")
    print("|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['scenario']} | {r['type']} | {r['fate']} | {r['basin_fraction']:.3f} | {r['active_nodes']} |")
