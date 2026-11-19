"""Custom exceptions for pypedigree."""

from __future__ import annotations

from typing import Any


class PedigreeError(Exception):
    """Base exception for all pypedigree errors."""


class PedigreeParseError(PedigreeError):
    """Raised when a pedigree file cannot be parsed."""


class PedigreeValidationError(PedigreeError):
    """Raised when pedigree data fails validation checks."""


class PedigreeLayoutError(PedigreeError):
    """Raised when the layout engine encounters an unresolvable situation."""


class PedigreeWarning(UserWarning):
    """Category for non-fatal pypedigree diagnostics.

    Warnings use this category so callers can filter them independently of
    other libraries' warnings, e.g.::

        warnings.simplefilter("error", PedigreeWarning)
    """


class PedigreeStyleWarning(PedigreeWarning):
    """Raised when a style definition contains something unusable."""


class PedigreeDataWarning(PedigreeWarning):
    """Raised when pedigree data is questionable but still renderable."""


class PedigreeBatchError(PedigreeError):
    """Raised by `render_families` when some families could not be drawn.

    The other families were drawn.

    Attributes:
        written: Family ID → path of each file written.
        failures: Family ID → error message, for each family not drawn.
    """

    def __init__(self, written: dict[str, Any], failures: dict[str, str]) -> None:
        self.written = written
        self.failures = failures
        listed = "; ".join(f"{fid}: {msg}" for fid, msg in failures.items())
        total = len(written) + len(failures)
        super().__init__(f"{len(failures)} of {total} families could not be drawn ({listed})")
