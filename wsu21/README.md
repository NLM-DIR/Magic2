# Exact matrix representations of sl(2) and sl(2|1)

A small set of Python programs that build the matrices of the Lie algebra sl(2)
and of the Lie superalgebra sl(2|1) in exact rational (and symbolic)
arithmetic, and verify every (super-)commutator, the grading operator, and the
Casimir and anticenter invariants.

Authors: Jean Thierry-Mieg (NLM/NIH) for the concepts
         Claude.AI:fable 5 for the programming and documentation
	 
License: public domain; no rights reserved.

## Contents

| File        | Role                                                              |
|-------------|-------------------------------------------------------------------|
| `matrix.py` | exact matrix class (library, not run directly)                    |
| `algebra.py`| a collection of generator matrices with parity (library)          |
| `sl2.py`    | sl(2) in the (a+1)-dimensional irreducible representation         |
| `sl21.py`   | sl(2\|1) Kac modules R(a,b) with numeric b, Casimirs, Matryoshka  |
| `sl21_b.py` | the same with b kept as a formal variable; zeta-Hermitian basis   |
| `MANUAL.md` | the mathematics and the command-line options                     |

All five Python files must sit in the same folder. There is nothing to compile and
nothing to install besides Python and SymPy.

## Requirements

- Python 3.9 or later
- SymPy (which pulls in its only dependency, mpmath, automatically)

Everything else used by the programs (`argparse`, `itertools`, `math`, ...)
is part of the Python standard library.

---

## Installation on Linux

### 1. Check that Python 3 is present

Almost every Linux distribution ships Python 3. In a terminal:

```bash
python3 --version
```

If this prints `Python 3.9` or higher, go to step 2. Otherwise install it with
your distribution's package manager:

```bash
# Debian, Ubuntu, Mint
sudo apt update
sudo apt install python3

# Fedora, RHEL, Rocky, Alma
sudo dnf install python3

# Arch, Manjaro
sudo pacman -S python
```

### 2. Install pip and the virtual-environment module

`pip` is Python's package installer. On Debian and Ubuntu it, and the `venv`
module, come as separate packages:

```bash
# Debian, Ubuntu, Mint
sudo apt install python3-pip python3-venv

# Fedora, RHEL, Rocky, Alma
sudo dnf install python3-pip

# Arch, Manjaro
sudo pacman -S python-pip
```

Check with `python3 -m pip --version`.

### 3. Install SymPy in a virtual environment (recommended)

Recent distributions (Ubuntu 23.04+, Debian 12+, Fedora 38+, Arch) refuse
`pip install` into the system Python and print
`error: externally-managed-environment`. The clean way around this is a
virtual environment: a private Python installation inside a folder, which
needs no administrator rights and cannot damage the system.

Go to the folder containing the `.py` files, then:

```bash
python3 -m venv .venv            # create the environment (once)
source .venv/bin/activate        # enter it (each new terminal)
pip install sympy                # install SymPy into it (once)
```

While the environment is active the prompt starts with `(.venv)`, and
`python` means the environment's Python. Leave it with `deactivate`.

### 3'. Alternatives, if you prefer not to use a virtual environment

Install SymPy from the distribution itself:

```bash
sudo apt install python3-sympy      # Debian, Ubuntu
sudo dnf install python3-sympy      # Fedora
sudo pacman -S python-sympy         # Arch
```

Distribution packages can lag a few versions behind, which is harmless here.
Or install for your user only (older distributions that still allow it):

```bash
python3 -m pip install --user sympy
```

### 4. Check the installation

```bash
python3 -c "import sympy; print(sympy.__version__)"
```

---

## Installation on macOS

macOS is Unix, and the Linux instructions apply almost unchanged. Install
Python 3 from <https://www.python.org/downloads/> or with Homebrew
(`brew install python`), then follow step 3 above (virtual environment and
`pip install sympy`).

---

## Installation on Windows

### 1. Install Python

Download the installer from <https://www.python.org/downloads/windows/> and run
it. **On the first screen, tick "Add python.exe to PATH"** before clicking
*Install Now*. pip is included automatically.

Alternatively, on Windows 10/11 you can install Python from the Microsoft Store,
or from a terminal with `winget install Python.Python.3.12`.

### 2. Install SymPy

Open *PowerShell* or *Command Prompt*, go to the folder containing the `.py`
files (`cd C:\path\to\folder`), and type:

```powershell
py -m pip install sympy
```

`py` is the Python launcher installed by the python.org installer; `python`
also works if you ticked "Add to PATH". A virtual environment works as on
Linux, with a different activation command:

```powershell
py -m venv .venv
.venv\Scripts\activate
pip install sympy
```

(If PowerShell refuses to run the activation script, use Command Prompt, or run
once `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`.)

### 3. Check the installation

```powershell
py -c "import sympy; print(sympy.__version__)"
```

### Note on special characters

The programs print a few non-ASCII symbols (✓, Σ, σ). They display correctly
in a modern Windows terminal. If you **redirect the output to a file**
(`> out.txt`) and Python stops with a `UnicodeEncodeError`, run it in UTF-8
mode:

```powershell
py -X utf8 sl21.py -a 1 -b 1 > out.txt
```

or set the variable once for the session: `set PYTHONUTF8=1` (Command Prompt)
or `$env:PYTHONUTF8=1` (PowerShell).

---

## First run

What the programs compute, and all their options, are described in
`MANUAL.md`.

In the folder containing the files (with the virtual environment active, if
you made one). Replace `python3` by `py` on Windows.

Each program prints its help when called without arguments:

```bash
python3 sl2.py
python3 sl21.py
python3 sl21_b.py
```

A few quick tests:

```bash
python3 sl2.py -a 2                  # sl(2) adjoint, 3x3
python3 sl21.py -a 0 -b 0            # sl(2|1) fundamental, 4x4
python3 sl21.py -a 1 -b 1            # sl(2|1) adjoint, 8x8
python3 sl21_b.py -a 0               # fundamental with formal b
python3 sl21_b.py -a 0 --zetaH --case 2
```

Each run prints the generator matrices and ends its verification sections with
lines such as

```
All sl(2|1) relations verified (nonzero and null). ✓
```

A line containing `FAILED` means a relation does not hold; that should never
happen with the distributed files.

### Running time

Plain runs take a second or two. The option `--casimirs` (Killing metric,
cubic tensor, Casimirs, anticenter) is much heavier: seconds for small `a` in
`sl21.py`, and noticeably longer in `sl21_b.py`, where b is symbolic. Large
`a` or Matryoshka depth `-N` increases the cost quickly; the programs refuse
sizes that would take unreasonably long.

## Keeping the output

The output is plain text; redirect it to a file to keep it:

```bash
python3 sl21.py -a 1 -b 1 --casimirs > R11.txt
```
