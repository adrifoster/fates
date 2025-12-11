"""
Builds/compiles any tests within the FATES repository
"""

import os
import shutil
from typing import Optional
from testing.framework.core.cime_interface import CIMEInterface

# # constants for this script
# _CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../cime")


class FortranTestBuilder:
    """Builder that prepares, configures, and compiles Fotran-based tests."""

    def __init__(
        self,
        *,
        build_dir: str,
        cmake_dir: str,
        mpilib: str = "mpi-serial",
        cime: CIMEInterface,
    ):
        self.build_dir = os.path.abspath(build_dir)
        self.cmake_dir = os.path.abspath(cmake_dir)
        self.mpilib = mpilib
        self.cime = cime

        self._pfunit_path: Optional[str] = None
        self._netcdf_c_path: Optional[str] = None
        self._netcdf_f_path: Optional[str] = None
        self._cmake_args: Optional[str] = None

    def build(self, *, make_j: int, clean: bool, verbose: bool):
        """Creates the build directory, resolves libraries, runs cmake,
        and compiles with make.

        Args:
            make_j (int): how many jobs to run in parallel for Make
            clean (bool): whether to clean the existing build before re-building
            verbose (bool): whether to use verbose printing in the compile commands
        """
        self._prepare_build_dir(clean)
        self._resolve_library_paths()
        self._run_cmake()
        self._run_make(make_j, clean, verbose)

    def _prepare_build_dir(self, clean: bool):
        """Preps build directory

        Args:
            clean (bool): whether or not to clean the build
        """
        if not os.path.isdir(self.build_dir):
            os.makedirs(self.build_dir)

        os.chdir(self.build_dir)

        if clean:
            self._clean_cmake_artifacts()

    def _clean_cmake_artifacts(self):
        """Deletes CMake-generated files in build directory."""
        if os.path.isfile("CMakeCache.txt"):
            os.remove("CMakeCache.txt")

        if os.path.isdir("CMakeFiles"):
            shutil.rmtree("CMakeFiles")

        for fname in os.listdir("."):
            if (
                fname in ("Macros.cmake", "env_mach_specific.xml")
                or fname.startswith("Depends")
                or fname.startswith(".env_mach_specific")
            ):
                try:
                    os.remove(fname)
                except IsADirectoryError:
                    shutil.rmtree(fname)

    def _resolve_library_paths(self):
        """Creates a fake CIME case to query PFUNIT/NETCDF paths and CMake args."""
        cime = self.cime

        machobj = cime.Machines()
        compiler = machobj.get_default_compiler()
        os_ = machobj.get_value("OS")

        # configure environment + macros
        cime.configure(
            machobj,
            self.build_dir,
            ["CMake"],
            compiler,
            self.mpilib,
            True,
            "nuopc",
            os_,
            unit_testing=True,
        )

        machspecific = cime.EnvMachSpecific(self.build_dir, unit_testing=True)
        fake_case = cime.FakeCase(compiler, self.mpilib, True, "nuopc", threading=False)
        machspecific.load_env(fake_case)

        args = [
            f"-DOS={os_}",
            f"-DMACH={machobj.get_machine_name()}",
            f"-DCOMPILER={compiler}",
            f"-DDEBUG={cime.stringify_bool(True)}",
            f"-DMPILIB={self.mpilib}",
            f"-Dcompile_threaded={cime.stringify_bool(False)}",
            f"-DCASEROOT={self.build_dir}",
        ]
        self._cmake_args = " ".join(args)

        self._pfunit_path = self._find_lib("PFUNIT_PATH")
        if "NETCDF" not in os.environ:
            self._netcdf_c_path = self._find_lib("NETCDF_C_PATH")
            self._netcdf_f_path = self._find_lib("NETCDF_FORTRAN_PATH")

    def _find_lib(self, varname: str):
        cime = self.cime
        tmp = cime.CmakeTmpBuildDir(macroloc=self.build_dir)
        with tmp:
            all_vars = tmp.get_makefile_vars(cmake_args=self._cmake_args).splitlines()
            for line in all_vars:
                if ":=" in line:
                    key, val = (s.strip() for s in line.split(":="))
                    if key == varname:
                        return val
        raise RuntimeError(f"{varname} not found during CIME library discovery.")

    def _run_cmake(self):
        """Run cmake"""
        cime = self.cime
        cmake_module_dir = os.path.join(cime.root, "CIME", "non_py", "src", "CMake")
        genf90_dir = os.path.join(cime.root, "CIME", "non_py", "externals", "genf90")

        cmd = [
            "cmake",
            "-C Macros.cmake",
            self.cmake_dir,
            f"-DCIMEROOT={cime.root}",
            f"-DSRC_ROOT={cime.get_src_root()}",
            f"-DCIME_CMAKE_MODULE_DIRECTORY={cmake_module_dir}",
            "-DCMAKE_BUILD_TYPE=CESM_DEBUG",
            f"-DCMAKE_PREFIX_PATH={self._pfunit_path}",
            "-DUSE_MPI_SERIAL=ON",
            "-DENABLE_GENF90=ON",
            f"-DCMAKE_PROGRAM_PATH={self._pfunit_path}",
        ]

        if self._netcdf_c_path:
            cmd.append(f"-DNETCDF_C_PATH={self._netcdf_c_path}")

        if self._netcdf_f_path:
            cmd.append(f"-DNETCDF_F_PATH={self._netcdf_f_path}")

        cmd.extend(self._cmake_args.split())
        cime.run_cmd_no_fail(" ".join(cmd), combine_output=True)

    def _run_make(self, make_j: int, clean: bool, verbose: bool):
        """Run make"""
        cime = self.cime

        if clean:
            cime.run_cmd_no_fail("make clean")

        cmd = ["make", "-j", str(make_j)]
        if verbose:
            cmd.append("VERBOSE=1")

        cime.run_cmd_no_fail(" ".join(cmd), combine_output=True)

    def build_exists(self, test_dir: str, test_exe: str | None = None) -> bool:
        """Checks to see if the build directory and any associated executables exist.

        Args:
        test_dir (str): test directory
        test_exe (str): test executable
        Returns:
        bool: whether or not build directory and associated executables exist
        """
        if not os.path.isdir(self.build_dir):
            return False

        if not os.path.isdir(os.path.join(self.build_dir, test_dir)):
            return False

        if test_exe and not os.path.isfile(
            os.path.join(self.build_dir, test_dir, test_exe)
        ):
            return False

        return True
