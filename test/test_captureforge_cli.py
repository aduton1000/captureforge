#!/usr/bin/env python3
"""test_captureforge_cli.py — tests for the unified `captureforge` CLI (stdlib-only).

Run:  python3 test/test_captureforge_cli.py

Checks: --help/--version + each subcommand --help work; `captureforge run --print-cmd`
constructs the expected `nextflow run main.nf …` command WITHOUT executing Nextflow;
the vendored-preferred + system-fallback resolution works BOTH with and without a system
`nextflow` on PATH; raw forms still importable.
"""
import os, subprocess, sys, tempfile, unittest
from pathlib import Path

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def _load_cli():
    if REPO not in sys.path:
        sys.path.insert(0, REPO)
    import captureforge_cli.cli as mod
    return mod

def cli(*args, env=None):
    e = dict(os.environ, PYTHONPATH=REPO)
    if env:
        e.update(env)
    p = subprocess.run([sys.executable, "-m", "captureforge_cli.cli", *args],
                       capture_output=True, text=True, cwd=REPO, env=e)
    return p.returncode, p.stdout, p.stderr


class CliHelp(unittest.TestCase):
    def test_top_help_lists_three_subcommands(self):
        rc, so, se = cli("--help")
        self.assertEqual(rc, 0)
        for c in ("run", "panel", "index"):
            self.assertIn(c, so)
        # no invented per-stage subcommands
        for bad in ("align", "anno", "burden", "cnv"):
            self.assertNotIn(f"\n  {bad} ", so)

    def test_version(self):
        rc, so, se = cli("--version")
        self.assertEqual(rc, 0); self.assertEqual(so.strip(), "0.1.0")

    def test_subcommand_helps(self):
        for c, tok in (("panel", "symbols"), ("run", "nextflow"), ("index", "fasta")):
            rc, so, se = cli(c, "--help")
            self.assertEqual(rc, 0, se)
            self.assertIn(tok, (so + se).lower())


class CliRunCmd(unittest.TestCase):
    def test_run_print_cmd_vendored_preferred(self):
        rc, so, se = cli("run", "--print-cmd", "-profile", "local,conda",
                         "-params-file", "params.yaml")
        self.assertEqual(rc, 0, se)
        self.assertIn("nextflow", so)
        self.assertIn("run", so)
        self.assertIn("main.nf", so)
        self.assertIn("-profile local,conda", so)
        # resolution is announced on stderr; vendored is preferred when present
        self.assertIn("nextflow resolved", se)
        self.assertIn(".tools/nextflow", se)

    def test_run_passthrough_and_pipeline_override(self):
        rc, so, se = cli("run", "--print-cmd", "--pipeline", "/tmp/x/main.nf",
                         "-profile", "test", "-resume")
        self.assertEqual(rc, 0, se)
        self.assertIn("/tmp/x/main.nf", so)
        self.assertIn("-resume", so)


class NextflowResolution(unittest.TestCase):
    """Exercise vendored-preferred + system-fallback both with and without system nextflow."""
    def setUp(self):
        self.cli = _load_cli()

    def test_vendored_preferred_even_with_system_present(self):
        with tempfile.TemporaryDirectory() as d:
            fake = Path(d) / "nextflow"; fake.write_text("#!/bin/sh\n"); fake.chmod(0o755)
            old = os.environ.get("PATH")
            os.environ["PATH"] = d + os.pathsep + (old or "")
            try:
                path, kind = self.cli._resolve_nextflow()
            finally:
                os.environ["PATH"] = old or ""
            # real repo has the vendored launcher -> it wins over the system one
            self.assertIn("vendored", kind)
            self.assertTrue(path.endswith(".tools/nextflow"))

    def test_falls_back_to_system_when_vendored_absent(self):
        with tempfile.TemporaryDirectory() as d:
            fake = Path(d) / "nextflow"; fake.write_text("#!/bin/sh\n"); fake.chmod(0o755)
            old_path = os.environ.get("PATH"); old_vendored = self.cli.VENDORED_NF
            os.environ["PATH"] = d + os.pathsep + (old_path or "")
            self.cli.VENDORED_NF = os.path.join(d, "no_such_nextflow")
            try:
                path, kind = self.cli._resolve_nextflow()
            finally:
                os.environ["PATH"] = old_path or ""; self.cli.VENDORED_NF = old_vendored
            self.assertIn("system", kind)
            self.assertEqual(path, str(fake))

    def test_none_when_neither_available(self):
        old_path = os.environ.get("PATH"); old_vendored = self.cli.VENDORED_NF
        os.environ["PATH"] = "/nonexistent_dir_for_test"
        self.cli.VENDORED_NF = "/no/such/nextflow"
        try:
            path, kind = self.cli._resolve_nextflow()
        finally:
            os.environ["PATH"] = old_path or ""; self.cli.VENDORED_NF = old_vendored
        self.assertIsNone(path); self.assertEqual(kind, "none")


class RawFormsStillWork(unittest.TestCase):
    def test_init_gene_list_importable_and_script_runs(self):
        # the wrapped helper is argv-callable AND the raw script still works
        mod_path = os.path.join(REPO, "bin", "init_gene_list.py")
        with tempfile.TemporaryDirectory() as d:
            genes = Path(d) / "g.txt"; genes.write_text("TP53\nEGFR\n")
            out = Path(d) / "panel.spec.csv"
            r = subprocess.run([sys.executable, mod_path, str(genes), "-o", str(out)],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertTrue(out.exists())
            self.assertIn("gene", out.read_text().splitlines()[0])


if __name__ == "__main__":
    unittest.main(verbosity=2)
