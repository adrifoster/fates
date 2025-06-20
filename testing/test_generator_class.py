"""Module for generating empty test classes"""
import os
import re
from abc import ABC, abstractmethod
from utils import snake_to_camel

_TEST_SUB_DIR = os.path.dirname(os.path.abspath(__file__))
_UNIT_TESTING_DIR = "unit_testing"
_FUNCTIONAL_TESTING_DIR = "functional_testing"
_TEMPLATE_DIR = os.path.join(_TEST_SUB_DIR, "templates")
_UNIT_CMAKE_TEMPLATE = "cmake_utest_template.txt"
_FUNCTIONAL_CMAKE_TEMPLATE = "cmake_ftest_template.txt"
_PFUNIT_TEMPLATE = "pfunit_template.txt"
_FUNCTIONAL_TEST_TEMPLATE = "fortran_test_template.txt"
_CLASS_TEMPLATE = "test_class_template.txt"
_TEST_CMAKELISTS = os.path.join(_TEST_SUB_DIR, "CMakeLists.txt")
_UNIT_TESTS_CONFIG = os.path.join(_TEST_SUB_DIR, "unit_tests.cfg")
_FUNCTIONAL_TESTS_CONFIG = os.path.join(_TEST_SUB_DIR, "functional_tests.cfg")
_LOAD_CLASS_FILE = os.path.join(_TEST_SUB_DIR, "load_functional_tests.py")


class GenerateTestClass(ABC):
    """Abstract base class for creating empty FATES tests"""

    def __init__(self, test_name: str, sub_dir: str = None):
        self.test_name = re.sub(r"_test$", "", test_name.lower())
        self.sub_dir = sub_dir
        self.test_dir = self.generate_test_dir()
        self.build_dir = self.generate_test_build_dir()
        self.module_name = snake_to_camel(self.test_name)
        self.file_name = self.generate_file_name()

        if not hasattr(self, "section_header"):
            raise ValueError("Subclass must define 'section_header'")
        
        if not hasattr(self, "config_file"):
            raise ValueError("Subclass must define 'config_file'")

    @abstractmethod
    def generate_test_dir(self) -> str:
        """Generates a test directory path (must be implemented by subclasses)."""

    @abstractmethod
    def generate_test_build_dir(self) -> str:
        """Generates a test build directory path (must be implemented by subclasses)."""

    @abstractmethod
    def generate_file_name(self) -> str:
        """ "Generates a file name for the test program (must be implemented by subclases)."""
        
    @abstractmethod
    def append_config_lines(self, lines):
        """Appends lines to the config file for the test (must be implemented by subclasses)."""

    @abstractmethod
    def create_cmake_file(self):
        """Creates the new CMakeLists.txt file for the test (must be implemented by subclasses)."""

    @abstractmethod
    def create_template_program(self):
        """Creates a template program file (must be implemented by subclasses)."""

    def create_test_directory(self):
        """Creates a new test directory

        Raises:
            FileExistsError: Directory already exists
        """
        if os.path.exists(self.test_dir):
            raise FileExistsError(f"Error: Directory '{self.test_dir}' already exists.")
        os.makedirs(self.test_dir)
        print(f"Created directory: {self.test_dir}")

    def update_main_cmake(self):
        """Updates the main testing/CMakeLists.txt file by adding a new test directory in
        the correct section.

        Raises:
            FileNotFoundError: testing/CMakeLists.txt file not found
        """
        if not os.path.exists(_TEST_CMAKELISTS):
            raise FileNotFoundError(f"Warning: {_TEST_CMAKELISTS} not found.")

        with open(_TEST_CMAKELISTS, "r", encoding="utf-8") as f:
            lines = f.readlines()

        new_entry = f"add_subdirectory({self.test_dir} {self.build_dir})\n"
        section_found = False
        updated_lines = []

        for line in lines:
            updated_lines.append(line)

            # insert new entry in the correct section
            if line.strip() == self.section_header:
                section_found = True

            # if we're in the correct section and hit another section header, insert it
            elif section_found and line.startswith("## "):
                updated_lines.insert(-1, new_entry)  # insert before new section header
                section_found = False  # stop inserting

        # if we never found a new section header, append at end
        if section_found:
            updated_lines.append(new_entry)

        # write modified contents back to file
        with open(_TEST_CMAKELISTS, "w", encoding="utf-8") as f:
            f.writelines(updated_lines)

        print(f"Updated {_TEST_CMAKELISTS} to include {self.test_dir}.")
        
    def update_config_file(self):
        """Loads config file lines"""
        
        if not os.path.exists(self.config_file):
            raise FileNotFoundError(f"ERROR: {self.config_file} not found.")
        
        with open(self.config_file, "r", encoding="utf-8") as f:
            lines = f.readlines()
            
        # strip trailing newlines and whitespace-only lines
        while lines and lines[-1].strip() == "":
            lines.pop()
        
        # add one blank line before new section
        lines.append("\n")
        
        # new config block
        lines = self.append_config_lines(lines)
        
        with open(self.config_file, "w", encoding="utf-8") as f:
            f.writelines(lines)

        print(f"Updated {self.config_file}.")        
        
    def setup_test(self):
        """Runs all setup steps: creates directory, updates CMake, and config."""

        self.create_test_directory()
        self.create_cmake_file()
        self.update_main_cmake()
        self.update_config_file()
        self.create_template_program()

    @staticmethod
    def load_template(template_name: str) -> str:
        """Load a template file from the templates directory

        Args:
            template_name (str): name of template

        Raises:
            FileNotFoundError: Template file not found

        Returns:
            str: template string
        """

        template_path = os.path.join(_TEMPLATE_DIR, template_name)
        if not os.path.exists(template_path):
            raise FileNotFoundError(f"ERROR: Template '{template_path}' not found.")
        with open(template_path, "r", encoding="utf-8") as f:
            return f.read()


class GenerateUnitTest(GenerateTestClass):
    """Class for generating empty unit tests"""

    section_header = "## Unit tests"
    config_file = _UNIT_TESTS_CONFIG

    def generate_test_dir(self) -> str:
        """Generates a test directory path

        Returns:
            str: test directory path
        """
        test_name = f"{self.test_name}_test"  # append "_test"
        if self.sub_dir:
            return os.path.join(_UNIT_TESTING_DIR, self.sub_dir, test_name)
        return os.path.join(_UNIT_TESTING_DIR, test_name)

    def generate_test_build_dir(self) -> str:
        """Generates a test build directory path

        Returns:
            str: test build directory path
        """
        return f"fates_{self.test_name}_utest"

    def generate_file_name(self) -> str:
        """Generates a file name for the test program.

        Returns:
            str: file name
        Raises:
            RuntimeError: self.module name not set yet
        """
        if not self.module_name:
            raise RuntimeError(
                "self.module_name must be set before file name generated!"
            )
        return f"test_{self.module_name}.pf"

    def create_cmake_file(self):
        """Creates the new CMakeLists.txt file for the test"""
        cmake_template = self.load_template(_UNIT_CMAKE_TEMPLATE)
        cmake_path = os.path.join(self.test_dir, "CMakeLists.txt")

        with open(cmake_path, "w", encoding="utf-8") as f:
            f.write(
                cmake_template.format(
                    file_name=self.file_name, module_name=self.module_name
                )
            )
        print(f"Created {cmake_path}.")

    def append_config_lines(self, lines):
        
        lines.append(f"[{self.test_name}]\n")
        lines.append(f"test_dir = {self.build_dir}\n")
        
        return lines
        
    def create_template_program(self):
        pfunit_template = self.load_template(_PFUNIT_TEMPLATE)
        pfunit_path = os.path.join(self.test_dir, self.file_name)
        
        with open(pfunit_path, "w", encoding="utf-8") as f:
            f.write(
                pfunit_template.format(
                    module_name=self.module_name
                )
            )
        print(f"Added template test files in {self.test_dir}.")


class GenerateFunctionalTest(GenerateTestClass):
    """Class for generating functional tests"""

    section_header = "## Functional tests"
    config_file = _FUNCTIONAL_TESTS_CONFIG

    def __init__(self, test_name, sub_dir=None):
        super().__init__(test_name, sub_dir)

        self.executable_name = f"{self.module_name}_exe"

    def generate_test_dir(self) -> str:
        """Generates a test directory path

        Returns:
            str: test directory path
        """
        if self.sub_dir:
            return os.path.join(_FUNCTIONAL_TESTING_DIR, self.sub_dir, self.test_name)
        return os.path.join(_FUNCTIONAL_TESTING_DIR, self.test_name)

    def generate_file_name(self) -> str:
        """Generates a file name for the test program.

        Returns:
            str: file name
        Raises:
            RuntimeError: self.module name not set yet
        """
        if not self.module_name:
            raise RuntimeError(
                "self.module_name must be set before file name generated!"
            )
        return f"Test{self.module_name}.F90"

    def generate_test_build_dir(self) -> str:
        """Generates a test build directory path

        Returns:
            str: test build directory path
        """
        return f"fates_{self.test_name}_ftest"

    def create_cmake_file(self):
        """Creates the new CMakeLists.txt file for the test"""
        cmake_template = self.load_template(_FUNCTIONAL_CMAKE_TEMPLATE)
        cmake_path = os.path.join(self.test_dir, "CMakeLists.txt")

        with open(cmake_path, "w", encoding="utf-8") as f:
            f.write(
                cmake_template.format(
                    file_name=self.file_name, executable_name=self.executable_name
                )
            )
        print(f"Created {cmake_path}.")
        
    def append_config_lines(self, lines):
        
        lines.append(f"[{self.test_name}]\n")
        lines.append(f"test_dir = {self.build_dir}\n")
        lines.append(f"test_exe = {self.executable_name}\n")
        lines.append("out_file = None\n")
        lines.append("use_param_file = False\n")
        lines.append("other_args = []\n")
        return lines

    def create_template_program(self):
        
        # create a fortran template file
        functional_test_template = self.load_template(_FUNCTIONAL_TEST_TEMPLATE)
        test_path = os.path.join(self.test_dir, self.file_name)
        
        with open(test_path, "w", encoding="utf-8") as f:
            f.write(
                functional_test_template.format(
                    module_name=self.module_name
                )
            )
        
        # create the python template file
        class_template = self.load_template(_CLASS_TEMPLATE)
        class_path = os.path.join(self.test_dir, f"{self.test_name}.py")
        
        with open(class_path, "w", encoding="utf-8") as f:
            f.write(
                class_template.format(
                    module_name=self.module_name
                )
            )
            
        print(f"Added template test files in {self.test_dir}.")
        
        # update the class load
        if not os.path.exists(_LOAD_CLASS_FILE):
            raise FileNotFoundError(f"ERROR: {_LOAD_CLASS_FILE} not found.")
        
        if self.sub_dir:
            name = f"{self.sub_dir}.{self.test_name}"
        else:
            name = self.test_name
            
        with open(_LOAD_CLASS_FILE, "a", encoding="utf-8") as f:
            f.write(f"from functional_testing.{name} import {self.module_name}\n")
        
        print(f"Updated {_LOAD_CLASS_FILE}.")
        

def generate_test(
    test_type: str, test_name: str, sub_dir: str = None
) -> GenerateTestClass:
    """Generates all boilerplate associated with a new test

    Args:
        test_type (str): test type ('unit' or 'functional')
        test_name (str): name of test
        sub_dir (str, optional): optional subdirectory to place. Defaults to None.

    Raises:
        RuntimeError: _description_

    Returns:
        GenerateTestClass: _description_
    """
    if test_type == "unit":
        return GenerateUnitTest(test_name, sub_dir)
    if test_type == "functional":
        return GenerateFunctionalTest(test_name, sub_dir)
    raise RuntimeError("test_type must be one of ['unit', 'functional']")
