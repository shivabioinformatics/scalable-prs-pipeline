#!/usr/bin/env python3
"""
store_results.py — Store PRS results in SQLite database

I'm using SQLite because it's a real SQL database in a single file.
No server setup, Python has it built in. In production this same schema
maps directly to PostgreSQL on Amazon RDS.
"""

import sys, os, sqlite3, argparse
from datetime import datetime

def parse_args():
    p = argparse.ArgumentParser(description="Store PRS results in SQLite")
    p.add_argument("--scores", required=True)
    p.add_argument("--gwas", required=True)
    p.add_argument("--qc-stats", required=True)
    p.add_argument("--db", required=True)
    p.add_argument("--run-id", default=None)
    return p.parse_args()

def create_schema(conn):
    """
    Three tables:
    - samples: patient scores + risk categories (what the clinician queries)
    - variants: GWAS weights (the lookup table — like DynamoDB in production)
    - qc_metrics: per-sample quality stats (flagging unreliable scores)
    """
    c = conn.cursor()
    c.execute("""CREATE TABLE IF NOT EXISTS samples (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        run_id TEXT NOT NULL, sample_id TEXT NOT NULL,
        prs_raw REAL NOT NULL, prs_zscore REAL NOT NULL,
        percentile REAL NOT NULL, risk_category TEXT NOT NULL,
        variants_used INTEGER, variants_imputed INTEGER,
        created_at TEXT NOT NULL, UNIQUE(run_id, sample_id))""")

    c.execute("""CREATE TABLE IF NOT EXISTS variants (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        rsid TEXT NOT NULL, chromosome TEXT, position INTEGER,
        effect_allele TEXT NOT NULL, other_allele TEXT,
        beta REAL NOT NULL, standard_error REAL,
        p_value REAL NOT NULL, allele_frequency REAL,
        UNIQUE(rsid))""")

    c.execute("""CREATE TABLE IF NOT EXISTS qc_metrics (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        run_id TEXT NOT NULL, sample_id TEXT NOT NULL,
        call_rate REAL, mean_depth REAL,
        variants_called INTEGER, variants_total INTEGER,
        pass_qc TEXT, UNIQUE(run_id, sample_id))""")

    # indexes for common queries
    c.execute("CREATE INDEX IF NOT EXISTS idx_risk ON samples(risk_category)")
    c.execute("CREATE INDEX IF NOT EXISTS idx_sid ON samples(sample_id)")
    c.execute("CREATE INDEX IF NOT EXISTS idx_rsid ON variants(rsid)")
    c.execute("CREATE INDEX IF NOT EXISTS idx_pval ON variants(p_value)")
    conn.commit()
    print("  Schema created: 3 tables + 4 indexes")

def load_scores(conn, path, run_id):
    c = conn.cursor()
    n = 0
    with open(path) as f:
        f.readline()
        for line in f:
            fs = line.strip().split("\t")
            if len(fs) < 7: continue
            c.execute("""INSERT OR REPLACE INTO samples
                (run_id,sample_id,prs_raw,prs_zscore,percentile,risk_category,variants_used,variants_imputed,created_at)
                VALUES (?,?,?,?,?,?,?,?,?)""",
                (run_id, fs[0], float(fs[1]), float(fs[2]), float(fs[3]), fs[4], int(fs[5]), int(fs[6]), datetime.now().isoformat()))
            n += 1
    conn.commit()
    print(f"  Inserted {n} scores into 'samples'")

def load_variants(conn, path):
    c = conn.cursor()
    n = 0
    with open(path) as f:
        f.readline()
        for line in f:
            fs = line.strip().split("\t")
            if len(fs) < 9: continue
            try:
                c.execute("""INSERT OR REPLACE INTO variants
                    (rsid,chromosome,position,effect_allele,other_allele,beta,standard_error,p_value,allele_frequency)
                    VALUES (?,?,?,?,?,?,?,?,?)""",
                    (fs[0], fs[1], int(fs[2]), fs[3], fs[4], float(fs[5]), float(fs[6]), float(fs[7]), float(fs[8])))
                n += 1
            except (ValueError, IndexError): continue
    conn.commit()
    print(f"  Inserted {n} variants into 'variants'")

def load_qc(conn, path, run_id):
    c = conn.cursor()
    n = 0
    in_section = False
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("## PER-SAMPLE QC"):
                in_section = True; continue
            if in_section and line.startswith("##"): break
            if in_section and line.startswith("sample\t"): continue
            if in_section and line:
                fs = line.split("\t")
                if len(fs) >= 6:
                    try:
                        c.execute("""INSERT OR REPLACE INTO qc_metrics
                            (run_id,sample_id,call_rate,mean_depth,variants_called,variants_total,pass_qc)
                            VALUES (?,?,?,?,?,?,?)""",
                            (run_id, fs[0], float(fs[1]), float(fs[2]), int(fs[3]), int(fs[4]), fs[5]))
                        n += 1
                    except (ValueError, IndexError): continue
    conn.commit()
    print(f"  Inserted {n} QC records into 'qc_metrics'")

def show_queries(conn):
    """example queries a clinician or researcher would run"""
    c = conn.cursor()
    print("\n  --- Example Queries ---")

    c.execute("SELECT risk_category, COUNT(*) FROM samples GROUP BY risk_category")
    print("\n  Risk distribution:")
    for r in c.fetchall(): print(f"    {r[0]}: {r[1]} patients")

    c.execute("SELECT sample_id, prs_zscore, percentile FROM samples ORDER BY prs_raw DESC LIMIT 5")
    print("\n  Top 5 highest risk:")
    for r in c.fetchall(): print(f"    {r[0]}: z={r[1]:.2f}, pctile={r[2]:.0f}")

    c.execute("SELECT rsid, chromosome, beta, p_value FROM variants ORDER BY p_value LIMIT 5")
    print("\n  Most significant variants:")
    for r in c.fetchall(): print(f"    {r[0]} (chr{r[1]}): beta={r[2]:.4f}, p={r[3]:.2e}")

def main():
    args = parse_args()
    run_id = args.run_id or datetime.now().strftime("RUN_%Y%m%d_%H%M%S")
    print("=" * 60)
    print("Storing PRS Results in Database")
    print(f"  Run ID: {run_id}")
    print(f"  Database: {args.db}")
    print("=" * 60)

    conn = sqlite3.connect(args.db)
    create_schema(conn)
    load_scores(conn, args.scores, run_id)
    load_variants(conn, args.gwas)
    load_qc(conn, args.qc_stats, run_id)
    show_queries(conn)
    conn.close()
    print(f"\nDone! Database: {args.db}")

if __name__ == "__main__":
    main()
