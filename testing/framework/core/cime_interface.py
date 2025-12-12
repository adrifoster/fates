"""Module for dealing with interfacing with CIME"""

import os
import sys


class CIMEInterface:
    """
    Centralized access to the CIME Python modules.

    This class:
    - Locates the CIME directory
    - Puts it on sys.path (once)
    - Imports the needed CIME modules
    - Exposes them as attributes

    Downstream code imports only this class, never CIME directly.
    """

    def __init__(self, fates_root: str):
        self.fates_root = os.path.normpath(fates_root)

        self._cime_path = None
        self.utils = None
        self.configure = None

        self._initialize()

    # -----------------------------
    # Internal helpers
    # -----------------------------

    def _initialize(self):
        """Locate CIME, set up sys.path, and import modules."""
        self._cime_path = self._find_cime()
        self._ensure_on_sys_path(self._cime_path)
        self._ensure_on_sys_path(os.path.join(self._cime_path, "CIME", "Tools"))

        # now that CIME paths are injected, import CIME modules
        self._import_cime_modules()

    def _find_cime(self) -> str:
        """Resolve path to cime/ directory relative to FATES root.

        Raises:
            RuntimeError: Cannot find CIME at expected location

        Returns:
            str: path to cime
        """
        candidate = os.path.normpath(os.path.join(self.fates_root, "../../cime"))
        if os.path.isdir(candidate):
            return candidate

        raise RuntimeError(f"Cannot find CIME at expected location: {candidate}")

    @staticmethod
    def _ensure_on_sys_path(path: str):
        """Insert a directory into sys.path only if not already present.

        Args:
            path (str): path to insert
        """
        if path not in sys.path:
            sys.path.insert(1, path)

    def _import_cime_modules(self):
        """Import the CIME modules needed by testing framework."""
        import CIME.utils
        import CIME.BuildTools.configure as cime_config
        import CIME.XML.machines as cime_machines
        import CIME.XML.env_mach_specific as env_mach_specific
        import CIME.build as cime_build

        self.configure = cime_config.configure
        self.FakeCase = cime_config.configure.FakeCase
        self.get_src_root = CIME.utils.get_src_root
        self.Machines = cime_machines.Machines
        self.EnvMachSpecific = env_mach_specific.EnvMachSpecific
        self.CmakeTmpBuildDir = cime_build.CmakeTmpBuildDir

    # -----------------------------
    # Public API
    # -----------------------------

    @property
    def root(self) -> str:
        """Public method to find CIME path

        Returns:
            str: path to CIME
        """
        return self._cime_path
