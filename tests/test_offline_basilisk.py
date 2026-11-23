import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class OfflineBasilisk(unittest.TestCase):
    def test_offline_default_build_targets(self):
        makefile = (ROOT / "basilisk/src/Makefile").read_text()
        targets = next(line.split() for line in makefile.splitlines()
                       if line.startswith("all: "))
        self.assertIn("qcc", targets)
        self.assertNotIn("bview2D", targets)
        self.assertNotIn("bview3D", targets)
        self.assertNotIn("bview2Dm", targets)

    def test_viewer_fails_closed_at_build_time(self):
        result = subprocess.run(
            [os.environ.get("CC", "cc"), "-E", "-x", "c",
             str(ROOT / "basilisk/src/display.h")],
            capture_output=True, check=False, text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Basilisk network display is disabled", result.stderr)
        self.assertNotIn("ws_socket_open", result.stdout)

    def test_installer_builds_local_tree_and_does_not_reset(self):
        with tempfile.TemporaryDirectory(prefix="Basilisk build ") as tmp:
            work = Path(tmp)
            installer = work / "reset_install_requirements.sh"
            shutil.copyfile(ROOT / installer.name, installer)
            src = work / "basilisk/src"
            src.mkdir(parents=True)
            marker = src / "keep-me"
            marker.touch()
            (src / "config").touch()
            (src / "Makefile").write_text("all:\n\t@touch built\n")
            shutil.copyfile(ROOT / "basilisk/src/display.h", src / "display.h")

            rejected = subprocess.run(
                ["bash", str(installer), "--hard"],
                cwd=work, capture_output=True, check=False,
            )
            self.assertNotEqual(rejected.returncode, 0)
            self.assertTrue(marker.exists())
            self.assertFalse((src / "built").exists())

            subprocess.run(["bash", str(installer)], cwd=work,
                           capture_output=True, check=True)
            self.assertTrue((src / "built").exists())
            self.assertTrue((work / ".project_config").is_file())
            self.assertTrue(marker.exists())

    def test_installer_does_not_build_unhardened_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            installer = work / "reset_install_requirements.sh"
            shutil.copyfile(ROOT / installer.name, installer)
            src = work / "basilisk/src"
            src.mkdir(parents=True)
            (src / "display.h").write_text("/* unmanaged vendor source */\n")
            (src / "Makefile").write_text("all:\n\t@touch built\n")
            rejected = subprocess.run(
                ["bash", str(installer)], cwd=work,
                capture_output=True, check=False, text=True,
            )
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn("hardened vendored Basilisk tree is missing",
                          rejected.stderr)
            self.assertFalse((src / "built").exists())


if __name__ == "__main__":
    unittest.main()
